function Q = findQ(data,freqMin,freqMax)
%index the frequencies and magnitudes
index = find(data(:,1)>freqMin & data(:,1)<freqMax);
freq = data(index,1);
magnitude = data(index,2);

%maximum freq
fmax = max(magnitude);
%range of the half max frequency
range = find(magnitude>max(magnitude)/2);
df = freq(range(end))-freq(range(1));

%Calculate Q
Q = fmax/df;