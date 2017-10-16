produce.two.dimensional.noise.local <-
    function( z
			,nr.frames = 60
            ,max.level = 1              # 0: global noise, >0: various degrees of localization
            ,win.type  = 'flat.hanning' # 'hanning' or 'flat.hanning'
            ,war.thr   = .1             
            ,overlap   = 40
            ,do.set.seed = FALSE
			,on.screen = FALSE
			) {
            
	#Produces a 2-dimensional correlated noise
	#Use the last observation as filter
    #Nested implementation to account for non-stationarities
	#Example:
	#z=Xln0[,,nr.obs]
	#At = produce.two.dimensional.noise.local(z=Xln0[,,nr.obs])
    #Created:
    #ned, September 2017
    
    ## Input image
    
    # make sure non-rainy pixels are set to zero
    min.value = min(z)
    orig.z = z
    z = z - min.value
    
    # store original field size
    orig.dim.x = dim(z)[1]
    orig.dim.y = dim(z)[2]
    orig.dm.size = c(orig.dim.x,orig.dim.y)
    
    # apply window to the image to limit spurious edge effects
    orig.window=build.two.dimensional.window(wsize=orig.dm.size,wtype=win.type)
    z = z*orig.window
    
    # now buffer the field with zeros to get a squared domain       <-- need this at the moment for the nested approach, but I guess we could try to avoid it
    dim.x = max(orig.dm.size) 
    dim.y = dim.x
    dm.size = c(dim.x,dim.y)
    ztmp = array(0,dim=dm.size)
    if(orig.dm.size[1] == dim.x){
        idx.buffer = round((dim.y - orig.dim.y)/2)
        ztmp[,idx.buffer:(idx.buffer+orig.dim.y-1)] = z
    }else{
        idx.buffer = round((dim.x - orig.dim.x)/2)
        ztmp[idx.buffer:(idx.buffer+orig.dim.x-1),] = z
    }
    z=ztmp
    
    ## Nested algorithm
       
    # prepare indices
    Idxi = matrix(c(1,dim.x),nrow=1) 
    Idxj = matrix(c(1,dim.y),nrow=1) 
    Idxipsd = matrix(c(1,2^max.level),nrow=1) 
    Idxjpsd = matrix(c(1,2^max.level),nrow=1) 
    
    # generate the FFT sample frequencies
    res.km = 1 
    freq = get.fftfreq(dim.x,res.km)
    fx = outer(freq*0, freq, FUN = '+')
    fy = outer(freq, freq*0, FUN = '+')
    freq.grid = sqrt(fx^2 + fy^2)
    
    # get global fourier filter
    mfilter0 = get.fourier.filter(z)
    # and allocate it to the final grid
    mfilter = array(0, dim=c(2^max.level,2^max.level,dim(mfilter0)[1],dim(mfilter0)[2]))
    mfilter = sweep(mfilter,c(3,4),mfilter0,'+')

    # now loop levels and build composite spectra
    level=0
    while(level < max.level){

        for(m in 1:dim(Idxi)[1]){
            # cat('---------\n')
            # print(c(level,m))
            # print(c(Idxi[m,],Idxj[m,]))
            # cat('------\n')
            
            # the indices of rainfall field
            Idxnext = split.field(Idxi[m,],Idxj[m,],2)
            Idxinext = Idxnext[[1]]; Idxjnext = Idxnext[[2]]
            # the indices of the field of fourier filters
            Idxpsdnext = split.field(Idxipsd[m,],Idxjpsd[m,],2)
            Idxipsdnext = Idxpsdnext[[1]]; Idxjpsdnext = Idxpsdnext[[2]]
            
            for(n in 1:dim(Idxinext)[1]){
                # print(c(n,Idxinext[n,],Idxjnext[n,]))

                mask = get.mask(dm.size[1],Idxinext[n,],Idxjnext[n,],win.type)
                war = sum((z*mask)>0)/(Idxinext[n,2]-Idxinext[n,1]+1)^2 
                if(war>war.thr){
                    # the new filter 
                    newfilter = get.fourier.filter(z*mask)
                    
                    # compute logistic function to define weights as function of frequency
                    # k controls the shape of the weighting function
                    merge.weights = logistic.function(1/freq.grid, k=0.05, x0 = (Idxinext[n,2] - Idxinext[n,1] + 1)/2)
                    newfilter = newfilter*(1-merge.weights)
                    
                    # perform the weighted average of previous and new fourier filters
                    if(length(dim(mfilter[Idxipsdnext[n,1]:Idxipsdnext[n,2],Idxjpsdnext[n,1]:Idxjpsdnext[n,2],,]))>2){
                        marg=c(3,4)
                    }else{
                        marg=c(1,2)
                    }
                    
                    mfilter[Idxipsdnext[n,1]:Idxipsdnext[n,2],Idxjpsdnext[n,1]:Idxjpsdnext[n,2],,]=sweep(mfilter[Idxipsdnext[n,1]:Idxipsdnext[n,2],Idxjpsdnext[n,1]:Idxjpsdnext[n,2],,],marg,merge.weights,'*')
                    
                    mfilter[Idxipsdnext[n,1]:Idxipsdnext[n,2],Idxjpsdnext[n,1]:Idxjpsdnext[n,2],,]=sweep(mfilter[Idxipsdnext[n,1]:Idxipsdnext[n,2],Idxjpsdnext[n,1]:Idxjpsdnext[n,2],,],marg,newfilter,'+')
                }
            }
        }
        
        # update indices
        level = level + 1
        Idx = split.field(c(1,dm.size[1]),c(1,dm.size[2]),2^level)
        Idxi = Idx[[1]]; Idxj = Idx[[2]]
        Idxpsd = split.field(c(1,2^max.level),c(1,2^max.level),2^level)
        Idxipsd  = Idxpsd [[1]]; Idxjpsd  = Idxpsd[[2]]

    }
    
    ## Power-filter images

	# produce normal noise array
    if(do.set.seed==TRUE) set.seed(42)
	white.noise = array(rnorm(dm.size[1]*dm.size[2]*nr.frames) 
                ,dim=c(dm.size[1] ,dm.size[2] ,nr.frames))
    
    # build composite image of correlated noise
    corr.noise = array(0,dim=c(dim.x,dim.y,nr.frames))
    sum.of.masks = array(0,dim=c(dim.x,dim.y,nr.frames))
    idxi = array(0,dim=c(2,1))
    idxj = array(0,dim=c(2,1))
    winsize = round( dm.size[1] / 2^max.level )
    
    # loop frames
	for(m in 1:nr.frames) {
        
        # get fourier spectrum of white noise field
        z.ftp = fft(white.noise[,,m])
        
        # loop rows
        for(i in 0:(2^max.level-1)){
            # loop columns
            for(j in 0:(2^max.level-1)){
            
            # apply fourier filtering with local filter
            this.filter = mfilter[i+1,j+1,,]
            z.ftp.fl = z.ftp * this.filter
            z.iftp = fft(z.ftp.fl ,inverse=T)
            z.cor = Re(z.iftp)
            
            # compute indices of local area
            idxi[1] = max( c(round(i*winsize - overlap/2 + 1), 1) )
            idxi[2] = min( c(round(idxi[1] + winsize + overlap/2 - 1), dm.size[1]) )
            idxj[1] = max( c(round(j*winsize - overlap/2 + 1), 1) )
            idxj[2] = min( c(round(idxj[1] + winsize + overlap/2 - 1), dm.size[2]) )
            
            # build mask and add local noise field to the composite image
            mask = get.mask(dm.size[1],idxi,idxj,win.type)
            corr.noise[,,m] = corr.noise[,,m] + z.cor*mask
            sum.of.masks[,,m] = sum.of.masks[,,m] + mask
            }
        }
	}	
    
    # normalize the sum
    idx = sum.of.masks > 0
    corr.noise[idx] = corr.noise[idx]/sum.of.masks[idx]
    
    # crop the image back to the original size
    difx = dim.x-orig.dim.x
    dify = dim.y-orig.dim.y
    output = corr.noise[(difx/2+1):(dim.x-difx/2),(dify/2+1):(dim.y-dify/2),]
    if(nr.frames==1) dim(output) = c(dim(output),1)

    # standardize the results to N(0,1)
    for(m in 1:nr.frames) {
        output[,,m]  = ( output[,,m] - mean(output[,,m]) ) / sd(output[,,m])
    }
    
    # plot results
    if(on.screen==TRUE) {
        x11()
        image(output[,,1], breaks=seq(-5,5,0.1), col = rainbow(100))
	}
	
	return(output)
}

build.two.dimensional.window <-
    function( wsize, wtype = 'flat.hanning'
            ) {
             
    two.dimensional.window <- tryCatch({
        # Asymmetric window
        # dim 1
        switch(wtype,
            hanning={
                w1d1 = 0.5 - 0.5*cos(2*pi*seq(wsize[1])/(wsize[1] - 1))},
            flat.hanning={
                T = wsize[1]/4
                W = wsize[1]/2
                B = seq(-W,W,length.out=2*W)
                R = abs(B)-T
                R[R<0]=0.
                A = 0.5*(1.0 + cos(pi*R/T))
                A[abs(B)>(2*T)]=0.0
                w1d1 = A},
            stop('Unknown window type.')
        ) 
        # dim 2
        switch(wtype,
            hanning={
                w1d2 = 0.5 - 0.5*cos(2*pi*seq(wsize[2])/(wsize[2] - 1))},
            flat.hanning={
                T = wsize[2]/4
                W = wsize[2]/2
                B = seq(-W,W,length.out=2*W)
                R = abs(B)-T
                R[R<0]=0.
                A = 0.5*(1.0 + cos(pi*R/T))
                A[abs(B)>(2*T)]=0.0
                w1d2 = A},
            stop('Unknown window type.')
        )
        
        # 2d window
        two.dimensional.window = sqrt(outer(w1d1,w1d2))
        },
        error=function(cond){
        # Symmetric window
        switch(wtype,
            hanning={
                w1d = 0.5 - 0.5*cos(2*pi*seq(wsize[1])/(wsize[1] - 1))},
            flat.hanning={
                T = wsize[1]/4
                W = wsize[1]/2
                B = seq(-W,W,length.out=2*W)
                R = abs(B)-T
                R[R<0]=0.
                A = 0.5*(1.0 + cos(pi*R/T))
                A[abs(B)>(2*T)]=0.0
                w1d = A},
            stop('Unknown window type.')
        )
        # 2d window
        two.dimensional.window = sqrt(outer(w1d,w1d))

        return(two.dimensional.window)
        }
    )  
    
    return(two.dimensional.window)
}

get.fourier.filter <-
    function( fieldin, do.norm = TRUE){
            
    # FFT of the field
    mfilter = fft(fieldin)
    
    # Normalize the real and imaginary parts
    if(do.norm==TRUE) {
        real.part = ( Re(mfilter) - mean(Re(mfilter)) ) / sd(Re(mfilter)) 
        imag.part = ( Im(mfilter) - mean(Im(mfilter)) ) / sd(Im(mfilter))
        mfilter = complex(real = real.part, imaginary = imag.part)
        mfilter = matrix(mfilter,nrow=dim(fieldin)[1],byrow=F)
    }
    
    # Extract the magnitude
	mfilter = Mod(mfilter)    
            
    return(mfilter)        
}
    
get.fftfreq <-
    function(n, d=1.0){
    
    if(n %% 2 == 0){
        f = c(0:(n/2-1), (-n/2):-1)/(d*n)
    }else{
        f = c(0:((n-1)/2), (-(n-1)/2):-1)/(d*n)
    }
    return(f)
}

split.field <-
    function(idxi,idxj,Segments){
    sizei = (idxi[2] - idxi[1]) + 1
    sizej = (idxj[2] - idxj[1]) + 1
    winsizei = round( sizei / Segments )
    winsizej = round( sizej / Segments )
    Idxi = array(0,dim=c(Segments^2,2))
    Idxj = array(0,dim=c(Segments^2,2))
    count=0
    for(i in seq(0,Segments-1)){
        for(j in seq(0,Segments-1)){
            count=count+1
            Idxi[count,1] = idxi[1] + i*winsizei
            Idxi[count,2] = min( c(Idxi[count,1] + winsizei - 1, idxi[2]) )
            Idxj[count,1] = idxj[1] + j*winsizej
            Idxj[count,2] = min( c(Idxj[count,1] + winsizej - 1, idxj[2]) )
        }
    }
    return(list(Idxi,Idxj))
}

get.mask <-
    function(Size,idxi,idxj,wintype){
    winsize = c(idxi[2] - idxi[1] + 1, idxj[2] - idxj[1] + 1)
    wind = build.two.dimensional.window(wsize = winsize, wtype = wintype)
    mask = array(0,dim=c(Size,Size)) 
    mask[idxi[1]:idxi[2],idxj[1]:idxj[2]] = wind
    return(mask)
}

logistic.function <-
    function(x,L = 1,k = 1,x0 = 0){
    return(L/(1 + exp(-k*(x - x0))))
}


