produce.two.dimensional.noise <-
   function( z
			,nr.frames = 60
            ,win.type = 'flat.hanning'
            ,do.set.seed = FALSE
			,on.screen = FALSE
			) {
	#Produces a 2-dimensional correlated noise
	#
	#Examples:
	#z=Xln0[,,nr.obs]
	#At = produce.two.dimensional.noise(z=Xln0[,,nr.obs])
    
    # make sure non-rainy pixels are set to zero
    min.value = min(z)
    orig.z = z
    z = z - min.value
    
    # store original field size
    orig.dim.x = dim(z)[1]
    orig.dim.y = dim(z)[2]
    orig.dm.size = c(orig.dim.x,orig.dim.y)
    
    # buffer the field with zeros to obtain a squared domain
    dim.x = max(orig.dm.size) 
    dim.y = dim.x
    dm.size = c(dim.x,dim.y)
    zs = array(0,dim=dm.size)
    if(orig.dm.size[1] == dim.x){
        idx.buffer = round((dim.y - orig.dim.y)/2)
        zs[,idx.buffer:(idx.buffer+orig.dim.y-1)] = z
    }else{
        idx.buffer = round((dim.x - orig.dim.x)/2)
        zs[idx.buffer:(idx.buffer+orig.dim.x-1),] = z
    }

	#Use the last observation as filter
    if(win.type == 'none'){
        mask = array(1,dim=dm.size)
    }else{
        mask = build.two.dimensional.window(wsize=dm.size,wtype=win.type)
    }
	mfilter = fft(zs*mask)
	mfilter = Mod(mfilter)
    	
	#Produce normal noise array
    if(do.set.seed==TRUE) set.seed(42)
	Ztf = array(rnorm(dm.size[1]*dm.size[2]*nr.frames) ,dim=c(dm.size[1],dm.size[2],nr.frames))
	
	#Power-filter images
	result = array(dim=c(orig.dim.x,orig.dim.y,nr.frames))
	for(m in 1:nr.frames) {
		z.ftp = fft(Ztf[,,m])

		#Multiply noise by the filter to build correlation
		z.ftp.fl = z.ftp * mfilter
		z.iftp = fft(z.ftp.fl ,inverse=T)
		#If everything is correct all imaginary parts should be 0
		z.cor = Re(z.iftp)
		
		#Crop a square equal to the original size
		#If smaller than that just stick the z.cor at the corner of a
		a = array(dim=c(orig.dim.x,orig.dim.y))
		if(orig.dim.y > dim(z.cor)[[1]]) {
			a[1:size,1:size] = z.cor
		} else {
			difx = dim.x-orig.dim.x
			dify = dim.y-orig.dim.y
			a = z.cor[(difx/2+1):(dim.x-difx/2),(dify/2+1):(dim.y-dify/2)]	
		}
		
		#Bring everything back to a N(0,1)
		a = (a-mean(a))/sd(a)
				
		if(on.screen==TRUE) {
			x11()
			image(a,breaks=seq(-5,5,0.1), col = rainbow(100))
		}
		result[,,m] = a
	}	
	
	
	return(result)
}

build.two.dimensional.window <-
    function( wsize, wtype = 'flat.hanning'
            ) {
             
    two.dimensional.window <- tryCatch(
        {
        # Asymmetric window
        
        # dim 1
        switch(wtype,
        
            hanning={
                w1d1 = 0.5 - 0.5*cos(2*pi*seq(wsize[1])/(wsize[1] - 1))
            },
            
            flat.hanning={
                T = wsize[1]/4
                W = wsize[1]/2
                B = seq(-W,W,length.out=2*W)
                R = abs(B)-T
                R[R<0]=0.
                A = 0.5*(1.0 + cos(pi*R/T))
                A[abs(B)>(2*T)]=0.0
                w1d1 = A
            },
            
            stop('Unsupported window type.')
        )
        
        # dim 2
        switch(wtype,
        
            hanning={
                w1d2 = 0.5 - 0.5*cos(2*pi*seq(wsize[2])/(wsize[2] - 1))
            },
            
            flat.hanning={
                T = wsize[2]/4
                W = wsize[2]/2
                B = seq(-W,W,length.out=2*W)
                R = abs(B)-T
                R[R<0]=0.
                A = 0.5*(1.0 + cos(pi*R/T))
                A[abs(B)>(2*T)]=0.0
                w1d2 = A
            },
            
            stop('Unsupported window type.')
        )
        
        # 2d window
        two.dimensional.window = sqrt(outer(w1d1,w1d2))
        
        },
        error=function(cond) {
        
        # Symmetric window
        switch(wtype,
        
            hanning={
                w1d = 0.5 - 0.5*cos(2*pi*seq(wsize[1])/(wsize[1] - 1))
            },
            
            flat.hanning={
                T = wsize[1]/4
                W = wsize[1]/2
                B = seq(-W,W,length.out=2*W)
                R = abs(B)-T
                R[R<0]=0.
                A = 0.5*(1.0 + cos(pi*R/T))
                A[abs(B)>(2*T)]=0.0
                w1d = A
            },
            
            stop('Unknown window type.')
        )
        # 2d window
        two.dimensional.window = sqrt(outer(w1d,w1d))

        return(two.dimensional.window)
        }
    )  
    
    return(two.dimensional.window)
}