function [filterA,filterB,center] = estfilt(nChannels,type,sRate,LowFreq,UpperFreq,dbg)

% This function returns the filter coefficients for a filter bank for a 
% given number of channels and with a variety of filter spacings. Also 
% returned are the filter centre frequencies.

% UpperFreq = 4999; LowFreq = 100;
% =========================================================================
% ------------------ Greenwood spacing of filters -------------------------
if strcmp(type,'greenwood') % ============= Greenwood spacing =============

	FS = sRate/2;
	nOrd = 6;
	[lower1,center,upper1] = greenwood(nChannels,LowFreq,UpperFreq,dbg);

	if FS<upper1(nChannels)
        useHigh = 1;
    else
        useHigh = 0;
	end

	filterA = zeros(nChannels,nOrd+1); 
    filterB = zeros(nChannels,nOrd+1);
    
    for i = 1:nChannels
        W1 = [lower1(i)/FS, upper1(i)/FS];
        if i == nChannels
            if useHigh == 0
                [b,a] = butter(3,W1);
            else
                [b,a] = butter(6,W1(1),'high');
            end
        else
            [b,a] = butter(3,W1);
        end
        filterB(i,1:nOrd+1) = b; % Save the coefficients 'b'.
        filterA(i,1:nOrd+1) = a; % Save the coefficients 'a'.
    end

% ---------------------- Linear filter spacing  ---------------------------
elseif strcmp(type,'linear') % ============== Linear spacing ==============

	FS = sRate/2;
	nOrd = 6;
	range = (UpperFreq-LowFreq);
	interval = range/nChannels;
	center = zeros(1,nChannels);

    % Figure out the center frequencies for all channels.
	for i = 1:nChannels  
		upper1(i) = LowFreq+(interval*i);
		lower1(i) = LowFreq+(interval*(i-1));
		center(i) = 0.5*(upper1(i)+lower1(i));
	end

	if FS < upper1(nChannels)
        useHigh = 1;
    else
        useHigh = 0;
	end

	filterA = zeros(nChannels,nOrd+1); 
	filterB = zeros(nChannels,nOrd+1); 

    for i = 1:nChannels
		W1 = [lower1(i)/FS, upper1(i)/FS];
		if i == nChannels
		  	if useHigh == 0
	    		[b,a] = butter(3,W1);
			else
				[b,a] = butter(6,W1(1),'high');
			end
	    else
		   	[b,a] = butter(3,W1);
	    end
	    filterB(i,1:nOrd+1) = b; % Save the coefficients 'b'.
	    filterA(i,1:nOrd+1) = a; % Save the coefficients 'a'.
    end

% ------------------- Logarithmic filter spacing  -------------------------
elseif strcmp(type,'log') % ================= Log spacing =================

	FS = sRate/2;
	nOrd = 6;
	range = log10(UpperFreq/LowFreq);
	interval = range/nChannels;
	center = zeros(1,nChannels);

    % Figure out the center frequencies for all channels.
	for i = 1:nChannels  
		upper1(i) = LowFreq*10^(interval*i);
		lower1(i) = LowFreq*10^(interval*(i-1));
		center(i) = 0.5*(upper1(i)+lower1(i));
	end

	if FS < upper1(nChannels)
        useHigh = 1;
    else
        useHigh = 0;
	end

	filterA = zeros(nChannels,nOrd+1); 
	filterB = zeros(nChannels,nOrd+1); 

	 for i = 1:nChannels
		W1 = [lower1(i)/FS, upper1(i)/FS];
		if i == nChannels
            if useHigh == 0
                  [b,a] = butter(3,W1);
            else
                 [b,a] = butter(6,W1(1),'high');
            end
	    else
	   	[b,a] = butter(3,W1);
	    end
	    filterB(i,1:nOrd+1) = b; % Save the coefficients 'b'.
	    filterA(i,1:nOrd+1) = a; % Save the coefficients 'a'.
    end

% =========================================================================
% --------------------- Shannon filter spacing  ---------------------------
elseif strcmp(type,'shannon') % ============ Shannon spacing ==============

    srat2 = sRate/2;
    rp = 1.5; % Passband ripple in dB.
    rs = 20;

% ----------- Pre-emphasis filter and Low-pass envelop filter -------------

    [bls,als] = ellip(1,rp,rs,1150/srat2,'high');
    [blo,alo] = butter(2,160/srat2);
  
    rs = 15.0;
    nOrd = 2; % Order of filter = 2*nOrd.
    nOrd2 = 2*nOrd+1; % Number of coefficients.
    nchan = nChannels;

    if nchan == 2
	
		filt2b = zeros(nchan,nOrd2);
		filt2a = zeros(nchan,nOrd2);
		[b,a] = ellip(nOrd,rp,rs,[50/srat2 1500/srat2]);
		filt2b(1,:) = b; 
        filt2a(1,:) = a;
		[b,a] = ellip(nOrd,rp,rs,[1500/srat2 4000/srat2]);
		filt2b(2,:) = b; 
        filt2a(2,:) = a;
	
        filtroA = zeros(nchan,nOrd2); 
        filtroB = zeros(nchan,nOrd2);
        filtroA = filt2a; 
        filtroB = filt2b;
        
    elseif nchan == 3
	
		filt3b = zeros(nchan,2*nOrd+1);
		filt3a = zeros(nchan,2*nOrd+1);
		crsf = [50 800 1500 4000];
        
		for i = 1:3
            lf = crsf(i)/srat2; 
            ef = crsf(i+1)/srat2;
            [b,a] = ellip(nOrd,rp,rs,[lf ef]);
            filt3b(i,:) = b; 
            filt3a(i,:) = a;
		end
	
        filtroA = zeros(nchan,2*nOrd+1); 
        filtroB = zeros(nchan,2*nOrd+1);
        filtroA = filt3a; 
        filtroB = filt3b;	 
        
    elseif nchan == 4
	
		filt4b = zeros(nchan,2*nOrd+1);
		filt4a = zeros(nchan,2*nOrd+1);
		crsf4 = [50 800 1500 2500 4000];
        
		for i = 1:4
            lf = crsf4(i)/srat2; 
            ef = crsf4(i+1)/srat2;
            [b,a] = ellip(nOrd,rp,rs,[lf ef]);
            filt4b(i,:) = b; 
            filt4a(i,:) = a;
		end
		
        filtroA = zeros(nchan,2*nOrd+1); 
        filtroB = zeros(nchan,2*nOrd+1);
        filtroA = filt4a; 
        filtroB = filt4b;
	
  end

% =========================================================================
% ------------------ Mel spacing of filters -------------------------------
elseif strcmp(type,'mel')  % ============= use Mel spacing ================

	FS = sRate/2;
	nOrd = 6;
	[lower1,center,upper1] = mel(nChannels,LowFreq,UpperFreq,dbg);

	if FS<upper1(nChannels)
        useHigh = 1;
    else
        useHigh = 0;
	end

    filterA = zeros(nChannels,nOrd+1);
    filterB = zeros(nChannels,nOrd+1);
    
    for i = 1:nChannels
        W1 = [lower1(i)/FS, upper1(i)/FS];
        if i == nChannels
            if useHigh == 0
                [b,a] = butter(3,W1);
            else
                [b,a] = butter(6,W1(1),'high');
            end
        else
            [b,a] = butter(3,W1);
        end
        filterB(i,1:nOrd+1) = b; % Save the coefficients 'b'.
        filterA(i,1:nOrd+1) = a; % Save the coefficients 'a'.
    end

else
    
    error('ERROR! filters must be log, greenwood, mel, linear or Shannon');

end
