produce.two.dimensional.noise <-
   function( z
			,size = 1024
			,nr.frames = 60
			,dim.x = 710
			,dim.y = 640
			,on.screen = FALSE
			) {
	#Produces a 2-dimensional correlated noise
	#
	#Examples:
	#z=Xln0[,,nr.obs]
	#At = produce.two.dimensional.noise(z=Xln0[,,nr.obs])

	#Use the last observation as filter
	dm.size = 512	
	zs = z[(floor((dim.x-dm.size)/2)+1):(floor((dim.x-dm.size)/2)+dm.size),(floor((dim.y-dm.size)/2)+1):((floor(dim.y-dm.size)/2)+dm.size)]	
	mfilter = fft(zs)
	mfilter = Mod(mfilter)
	tmp = array(0,dim=c(1024,1024))
	tmp[1:256,1:256] = mfilter[1:256,1:256]
	tmp[(1024-256):1024,1:256] = mfilter[256:512,1:256]
	tmp[1:256,(1024-256):1024] = mfilter[1:256,256:512]
	tmp[(1024-256):1024,(1024-256):1024] = mfilter[256:512,256:512]
	mfilter = tmp
	
	#Produce normal noise array
	Ztf = array(rnorm(size^2*nr.frames) ,dim=c(size ,size ,nr.frames))
	
	#Power-filter images
	result = array(dim=c(710,640,nr.frames))
	for(m in 1:nr.frames) {
		z.ftp = fft(Ztf[,,m])

		#Multiply noise by the filter to build correlation
		z.ftp.fl = z.ftp * mfilter
		z.iftp = fft(z.ftp.fl ,inverse=T)
		#If everything is correct all imaginary parts should be 0
		z.cor = Re(z.iftp)
		
		#Crop a square equal to 710x640
		#If smaller than that just stick the z.cor at the corner of a
		a = array(dim=c(dim.x,dim.y))
		if(dim.y > dim(z.cor)[[1]]) {
			a[1:size,1:size] = z.cor
		} else {
			difx = size-dim.x
			dify = size-dim.y
			a = z.cor[(difx/2+1):(dim.x+difx/2),(dify/2+1):(dim.y+dify/2)]	
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
