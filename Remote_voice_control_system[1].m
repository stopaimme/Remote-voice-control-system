%% 信号采集
clc
clear all
close all
recObj = audiorecorder;
Fs = 16000 ; 
nBits = 8 ; 
nChannels = 1 ; 
ID = -1; % default audio input device 
recObj = audiorecorder(Fs,nBits,nChannels,ID);%采集到的信号
disp('Start speaking.')
recordblocking(recObj,1);
disp('End of Recording.');
myRecording = getaudiodata(recObj);
figure(1);
subplot(3,1,1);
plot(myRecording);
title('原始声音信号')
xlabel('times')
ylabel('amplitude')

%% 低通滤波
FIR=load('FIR.mat');%调用滤波器
lowpassRecording=filter(FIR.Num,1,myRecording);%低通滤波
N=length(myRecording); %长度
n=0:N-1;
Y1=fft(myRecording); %对原始信号做FFT变换
Y2=fft(lowpassRecording);%对滤波信号做FFT变换
P2 = abs(Y1/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);

P2_L = abs(Y2/N);
P1_L = P2_L(1:N/2+1);
P1_L(2:end-1) = 2*P1_L(2:end-1);
f = Fs*(0:(N/2))/N;
subplot(3,1,2);
plot(f,P1) 
title('原始声音信号频谱')
xlabel('frequence(Hz)')
ylabel('amplitude')
subplot(3,1,3);
plot(f,P1_L) 
title('滤波后声音信号频谱')
xlabel('frequence(Hz)')
ylabel('amplitude')

%% 信号量化
partition = [-1:1/128:1];     %量化范围与间隔
codebook = [-1-1/128:1/128:1]; 
[index,out]= quantiz(lowpassRecording,partition,codebook);  %信号量化
figure(2);
plot(n/N,lowpassRecording);
hold on
plot(n/N,out);
hold off
title('Quantization of origin signal ')
xlabel('Time')
ylabel('Amplitude')
legend('Original sampled signal','Quantized signal');

%% 信源编码
a = zeros(1,256);
for i= 1:16000 %统计每个符号出现次数
    a(index(i))=a(index(i))+1;
end
symbols = 0:255; 
    a = a/16000;   %计算每个符号概率
dict = huffmandict(symbols,a);   %生成码书
code = huffmanenco(index,dict);   %进行编码
disp('原始信号长度')
disp(length(index))
disp('编码信号长度')
disp(length(code))
disp('压缩率')
disp(8*length(index)/length(code))

%% 信道编码(7 4循环码)
len = length(code);    
rec = mod(len,4);
newlen = ceil(len/4)*4;  %对原始信号补0，以实现被4整除
newcode=zeros(1,newlen);
newcode(1:len)=code;
newcode=reshape(newcode,[],4);
disp('生成多项式')
pol = cyclpoly(7,4)  %获得生成多项式
[parmat,genmat] = cyclgen(7,pol);%获得生成矩阵和监督矩阵
disp('生成矩阵')
genmat
disp('监督矩阵')
parmat
channelenc = encode(newcode,7,4,'linear/binary',genmat); %信道编码
channelenc=reshape(channelenc,1,[]);
clen = length(channelenc);
re = mod(clen,2);
cllen = ceil(clen/2)*2;
nnewcode=zeros(1,cllen);
nnewcode(1:clen)=channelenc;

nnewcode = reshape(nnewcode,[],1);
cllen = cllen/2;
qpskcode = zeros(1,cllen);  %对信道编码的结果进行2bit分组，以实现QPSK
    for i=1:cllen
        qpskcode(i)=nnewcode(2*i-1)+nnewcode(2*i)*2;
    end

%% 信号调制 qpsk
M = 4; %modulation order

Txmsg = pskmod(qpskcode,M,pi/M);  %获得发送信号
Txmsg = awgn(Txmsg,10);      %信号加噪
scatterplot(Txmsg)
%% 信号解调 qpsk
Rxmsg = pskdemod(Txmsg,4,pi/4);   %发送信号解调
lenrx = length(Rxmsg);
iqpsk = zeros(1,2*lenrx);
for i=1:lenrx               %把4进制信号恢复成2进制信号
        iqpsk(2*i)=floor(Rxmsg(i)/2);
        iqpsk(2*i-1)=mod(Rxmsg(i),2);
end
reshape(iqpsk,1,[]);
decqpsk=iqpsk(1:2*lenrx-re);   %得到原始的信道编码结果
disp('信道传输误码率')
[number,ratio] = biterr(channelenc,decqpsk)
%% 信道解调
decqpsk = reshape(decqpsk,[],7);
channeldec = decode(decqpsk,7,4,'linear/binary',genmat); %信道译码
newcode=reshape(newcode,[],1);
channeldec=reshape(channeldec,[],1);
%% 信源解调
if(mod(rec,4)~=0)
sourceenc = channeldec(1:length(channeldec)-(4-mod(rec,4)));  %得到原始的信源编码结果
else
   sourceenc = channeldec(1:length(channeldec));    
end
%% 信号重建
out = huffmandeco(sourceenc, dict);    %信源译码
disp('Huffman编码误码率')
[number,ratio] = biterr(code,sourceenc)
out = (out-128)/128;   %信号重建
%% 重建后信号和原始信号比较
figure(4);
subplot(2,1,1);
plot(out) 
xlim([0,16000]);
title('重建后信号')
xlabel('frequence(Hz)')
ylabel('amplitude')
subplot(2,1,2);
plot(myRecording) 
title('原始信号')
xlabel('frequence(Hz)')
ylabel('amplitude')

Y3=fft(out);
P3 = abs(Y3/N);
P3_L = P3(1:N/2+1);
P3_L(2:end-1) = 2*P3_L(2:end-1);
figure(5);
subplot(2,1,1);
plot(P1) 
xlim([0,8000]);
title('原始声音信号频谱')
xlabel('frequence(Hz)')
ylabel('amplitude')
subplot(2,1,2);
plot(P3) 
xlim([0,8000]);
title('重建后声音信号频谱')
xlabel('frequence(Hz)')
ylabel('amplitude')





