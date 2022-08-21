DFs = 10000;
[y,Fs] = audioread('Sound.wav');
y = (y(:,1) + y(:,2)) / 2;
[p, q] = rat(DFs / Fs);
resampled_y = resample(y, p, q);
range = max(abs(min(resampled_y)), max(resampled_y) * 128 / 127);
modified_y = resampled_y / range;
quantized_y = cast(modified_y * 128 + 128, 'uint8');
binary_y = reshape(de2bi(quantized_y)',1,[]);
fc1 = 80000;
fc2 = 160000;
df = fc2 - fc1;
fs = 400000;
output = zeros(1,length(binary_y) * fs / df);
time = 0:1 / fs:1 / df - 1 / fs;
phi = rand() * 2 * pi;
for i = 1:length(binary_y)
    if binary_y(i) == 0
        output((i - 1) * fs / df + 1:i * fs / df) = sqrt(2) * cos(2 * pi * fc1 * time + phi);
    else
        output((i - 1) * fs / df + 1:i * fs / df) = sqrt(2) * cos(2 * pi * fc2 * time + phi);
    end
end
c1 = sqrt(2) * cos(2 * pi * fc1 * time);
s1 = sqrt(2) * sin(2 * pi * fc1 * time);
c2 = sqrt(2) * cos(2 * pi * fc2 * time);
s2 = sqrt(2) * sin(2 * pi * fc2 * time);
SNRs = -15:5:15;
received = zeros(length(SNRs),length(output));
binary_r = zeros(length(SNRs),length(binary_y));
for snr = SNRs
    i = (snr - SNRs(1)) / 5 + 1;
    received(i,:) = awgn(output, snr);
    for j = 1:length(binary_y)
        rc1 = dot(received(i,(j - 1) * fs / df + 1:j * fs / df), c1);
        rs1 = dot(received(i,(j - 1) * fs / df + 1:j * fs / df), s1);
        rc2 = dot(received(i,(j - 1) * fs / df + 1:j * fs / df), c2);
        rs2 = dot(received(i,(j - 1) * fs / df + 1:j * fs / df), s2);
        if rc1 ^ 2 + rs1 ^ 2 >= rc2 ^ 2 + rs2 ^ 2
            binary_r(i,j) = 0;
        else
            binary_r(i,j) = 1;
        end
    end
    reconstructed_y = (cast(bi2de(reshape(binary_r(i,:),8,[])'), 'double') - 128) / 128;
    audiowrite("reconstructed_y with SNR = " + snr + ".wav", reconstructed_y, DFs);
end
binary_r = cast(binary_r, 'uint8');
pr = zeros(1,length(SNRs));
for i = 1:length(pr)
   pr(i) = sum(abs(binary_y - binary_r(i,:))) / length(binary_y);
end
%%
DFs = 10000;
[y,Fs] = audioread('Sound.wav');
y = (y(:,1) + y(:,2)) / 2;
[p, q] = rat(DFs / Fs);
resampled_y = resample(y, p, q);
range = max(abs(min(resampled_y)), max(resampled_y) * 32 / 31);
modified_y = resampled_y / range;
byte_y = de2bi(cast(modified_y * 32 + 32, 'uint8'));
additional = ceil(length(byte_y) / 33) * 33 - length(byte_y);
quantized_y = [byte_y; repmat([0 0 0 0 0 1],additional,1)];
binary_y = reshape(quantized_y',1,[]);
gf_y = gf(reshape(quantized_y',99,[])');
bch_y = bchenc(gf_y,127,99);
binary_bch_y = reshape(bch_y.x',1,[]);
fc1 = 80000;
fc2 = 160000;
df = fc2 - fc1;
fs = 400000;
output = zeros(1,length(binary_bch_y) * fs / df);
time = 0:1 / fs:1 / df - 1 / fs;
phi = rand() * 2 * pi;
for i = 1:length(binary_bch_y)
    if binary_bch_y(i) == 0
        output((i - 1) * fs / df + 1:i * fs / df) = sqrt(2) * cos(2 * pi * fc1 * time + phi);
    else
        output((i - 1) * fs / df + 1:i * fs / df) = sqrt(2) * cos(2 * pi * fc2 * time + phi);
    end
end
c1 = sqrt(2) * cos(2 * pi * fc1 * time);
s1 = sqrt(2) * sin(2 * pi * fc1 * time);
c2 = sqrt(2) * cos(2 * pi * fc2 * time);
s2 = sqrt(2) * sin(2 * pi * fc2 * time);
SNRs = 0:5:15;
received = zeros(length(SNRs),length(output));
binary_bch_r = zeros(length(SNRs),length(binary_bch_y));
binary_r = zeros(length(SNRs),length(binary_y));
for snr = SNRs
    i = (snr - SNRs(1)) / 5 + 1;
    received(i,:) = awgn(output, snr);
    for j = 1:length(binary_bch_y)
        rc1 = dot(received(i,(j - 1) * fs / df + 1:j * fs / df), c1);
        rs1 = dot(received(i,(j - 1) * fs / df + 1:j * fs / df), s1);
        rc2 = dot(received(i,(j - 1) * fs / df + 1:j * fs / df), c2);
        rs2 = dot(received(i,(j - 1) * fs / df + 1:j * fs / df), s2);
        if rc1 ^ 2 + rs1 ^ 2 >= rc2 ^ 2 + rs2 ^ 2
            binary_bch_r(i,j) = 0;
        else
            binary_bch_r(i,j) = 1;
        end
    end
    binary_r(i,:) = reshape(bchdec(gf(reshape(binary_bch_r(i,:),127,[])'),127,99).x',1,[]);
    reconstructed_y = (cast(bi2de(reshape(binary_r(i,:),6,[])'), 'double') - 32) / 32;
    audiowrite("(BCH) reconstructed_y with SNR = " + snr + ".wav", reconstructed_y, DFs);
end
binary_r = cast(binary_r, 'uint8');
pr = zeros(1,length(SNRs));
for i = 1:length(pr)
   pr(i) = sum(abs(binary_y - binary_r(i,:))) / length(binary_y);
end
%%
DFs = 10000;
[y,Fs] = audioread('Sound.wav');
y = (y(:,1) + y(:,2)) / 2;
[p, q] = rat(DFs / Fs);
resampled_y = resample(y, p, q);
MUs = [10, 127, 255];
for mu = MUs
    range = max(abs(min(resampled_y)), max(resampled_y) / compand(31 / 32,mu,1,'mu/expander'));
    modified_y = compand(resampled_y / range,mu,1);
    byte_y = de2bi(cast(modified_y * 32 + 32, 'uint8'));
    additional = ceil(length(byte_y) / 33) * 33 - length(byte_y);
    quantized_y = [byte_y; repmat([0 0 0 0 0 1],additional,1)];
    binary_y = reshape(quantized_y',1,[]);
    gf_y = gf(reshape(quantized_y',99,[])');
    bch_y = bchenc(gf_y,127,99);
    binary_bch_y = reshape(bch_y.x',1,[]);
    fc1 = 80000;
    fc2 = 160000;
    df = fc2 - fc1;
    fs = 400000;
    output = zeros(1,length(binary_bch_y) * fs / df);
    time = 0:1 / fs:1 / df - 1 / fs;
    phi = rand() * 2 * pi;
    for i = 1:length(binary_bch_y)
        if binary_bch_y(i) == 0
            output((i - 1) * fs / df + 1:i * fs / df) = sqrt(2) * cos(2 * pi * fc1 * time + phi);
        else
            output((i - 1) * fs / df + 1:i * fs / df) = sqrt(2) * cos(2 * pi * fc2 * time + phi);
        end
    end
    c1 = sqrt(2) * cos(2 * pi * fc1 * time);
    s1 = sqrt(2) * sin(2 * pi * fc1 * time);
    c2 = sqrt(2) * cos(2 * pi * fc2 * time);
    s2 = sqrt(2) * sin(2 * pi * fc2 * time);
    SNRs = 0:5:15;
    received = zeros(length(SNRs),length(output));
    binary_bch_r = zeros(length(SNRs),length(binary_bch_y));
    binary_r = zeros(length(SNRs),length(binary_y));
    for snr = SNRs
        i = (snr - SNRs(1)) / 5 + 1;
        received(i,:) = awgn(output, snr);
        for j = 1:length(binary_bch_y)
            rc1 = dot(received(i,(j - 1) * fs / df + 1:j * fs / df), c1);
            rs1 = dot(received(i,(j - 1) * fs / df + 1:j * fs / df), s1);
            rc2 = dot(received(i,(j - 1) * fs / df + 1:j * fs / df), c2);
            rs2 = dot(received(i,(j - 1) * fs / df + 1:j * fs / df), s2);
            if rc1 ^ 2 + rs1 ^ 2 >= rc2 ^ 2 + rs2 ^ 2
                binary_bch_r(i,j) = 0;
            else
                binary_bch_r(i,j) = 1;
            end
        end
        binary_r(i,:) = reshape(bchdec(gf(reshape(binary_bch_r(i,:),127,[])'),127,99).x',1,[]);
        reconstructed_y = compand((cast(bi2de(reshape(binary_r(i,:),6,[])'), 'double') - 32) / 32,mu,1,'mu/expander');
        audiowrite("(mu_law = " + mu + ") reconstructed_y with SNR = " + snr + ".wav", reconstructed_y, DFs);
    end
end