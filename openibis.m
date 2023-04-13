% Original script from the work of:
% Connor, C. W. (2022). Open Reimplementation of the BIS Algorithms for Depth of Anesthesia. Anesthesia & Analgesia, 135(4), 855-864.
% It has been slightly modified since the original publication had a few 'end' statements missing.


% Main function
function depthOfAnesthesia = openibis(eeg)
    [Fs, stride]    =   deal(128,0.5);
    [BSRmap, BSR]   = suppression(eeg, Fs, stride);
    components      = logPowerRatios(eeg, Fs, stride, BSRmap);
    depthOfAnesthesia = mixer(components, BSR);
end

function [BSRmap, BSR]  = suppression(eeg, Fs, stride)
    [N, nStride]    = nEpochs(eeg, Fs, stride);
    BSRmap          = zeros(N, 1);
    for n = 1:N
        x = segment(eeg, n + 6.5, 2, nStride);
        BSRmap(n) = all(abs(x-baseline(x)) <= 5);
    end
    BSR         = 100 * movmean(BSRmap, [(63/stride) - 1, 0]);
end

function components = logPowerRatios(eeg, Fs, stride, BSRmap)
    [N, nStride]      = nEpochs(eeg, Fs, stride);
    [B, A]            = butter(2, 0.65/(Fs/2), 'high');
    eegHiPassFiltered = filter(B, A, eeg);
    psd               = nan(N, 4 * nStride/2);
    suppressionFilter = piecewise(0:0.5:63.5, [0, 3, 6], [0, 0.25, 1]).^2;
    components        = nan(N, 3); 
    for n = 1:N
        if isNotBurstSuppressed(BSRmap, n, 4)
            psd(n,:) = powerSpectralDensity(segment(eegHiPassFiltered, n + 4, 4, nStride)); 
            if sawtoothDetector(segment(eeg, n+4, 4, nStride), nStride)
                psd(n, :) = suppressionFilter .* psd(n, :); 
            end
        end
        
        thirtySec   = timeRange(30, n, stride); 
        vhighPowerConc = sqrt(mean(psd(thirtySec, bandRange(39.5, 46.5, 0.5)) .* psd(thirtySec, bandRange(40, 47, 0.5)), 2 )); 
        wholePowerConc = sqrt(mean(psd(thirtySec, bandRange(0.5, 46.5, 0.5)) .* psd(thirtySec, bandRange(1, 47, 0.5)), 2 )); 
        midBandPower    = prctmean(nanmean(10*log10(psd(thirtySec, bandRange(11, 20, 0.5))),1),50, 100);
        components(n ,1) = meanBandPower(psd(thirtySec, :), 30, 47, 0.5) - midBandPower; 
        components(n, 2) = trimmean (10 * log10(vhighPowerConc./ wholePowerConc), 50);
        components(n, 3) = meanBandPower(psd(thirtySec, :), 0.5, 4, 0.5) - midBandPower;
        
    end
end

function y = powerSpectralDensity(x)
    f       = fft(blackman(length(x)) .* (x - baseline(x)));
    y       = 2 * abs(f(1:length(x)/2)').^2 / (length(x)*sum(blackman(length(x)).^2));
end


function y = sawtoothDetector(eeg, nStride)
    saw     = [zeros(1, nStride -5) 1:5]'; 
    saw     = (saw - mean(saw)) / std(saw, 1); 
    r       = 1:(length(eeg) - length(saw)); 
    v       = movvar(eeg, [0 length(saw)-1], 1);
    m       = ([conv(eeg, flipud(saw), 'valid') conv(eeg, saw, 'valid')] / length(saw)) .^2;
    y       = max([(v(r)>10).*m(r, 1)./v(r); (v(r)>10).*m(r, 2) ./v(r) ]) >0.63;
end

function y = mixer(components, BSR)
    sedationScore   = scurve(components(:, 1), 104.4, 49.4, -13.9, 5.29);
    generalScore    = piecewise(components(:,2), [-60.89, -30], [-40, 43-1]);
    generalScore    = generalScore + scurve(components(:, 2), 61.3, 72.6, -24.0, 3.55) .*(components(:,2)>= -30);
    bsrScore        = piecewise(BSR, [0, 100], [50, 0]);
    generalWeight   = piecewise( components(:,3), [0,5], [0.5, 1]) .* (generalScore < sedationScore);
    bsrWeight       = piecewise(BSR, [10, 50], [0, 1]);
    x               = (sedationScore .* (1-generalWeight)) + (generalScore .* generalWeight);
    y               = piecewise(x, [-40, 10, 97, 110], [0, 10, 97, 100]).*(1-bsrWeight) + bsrScore.*bsrWeight; 
end

% Helper functions

function [N, nStride]    = nEpochs(eeg, Fs, stride), nStride = Fs*stride; N=floor((length(eeg )-Fs )/ nStride)-10;end
function y = meanBandPower(psd, from, to, bins), v = psd(:, bandRange(from,to,bins)); y=mean(10*log10(v(~isnan(v)))); end 
function y = bandRange(from, to, bins), y = ((from/bins):(to/bins)) +1; end
function y = baseline(x),   v = (1:length(x))' .^(0:1); y = v * (v\x); end
function y = bound(x, lowerBound, upperBound), y = min(max(x, lowerBound), upperBound);end
function y = segment(eeg, from, number, nStride), y = eeg(from*nStride + (1:number*nStride)); end
function y = isNotBurstSuppressed(BSRmap, n, p), y = ~((n<p) || any(BSRmap(n+((1-p):0)))); end
function y = timeRange(seconds, n, stride),      y = max(1, (n -(seconds/stride)+1)):n;end
function y = prctmean(x, lo, hi),   v = prctile(x, [lo hi]);  y = mean(x(x>=v(1) & x<=v(2))); end
function y = piecewise(x, xp, yp),   y = interp1(xp, yp, bound(x, xp(1), xp(end))); end
function y = scurve(x, Eo, Emax, x50, xwidth), y = Eo - Emax./(1+exp((x-x50)/xwidth)); end


