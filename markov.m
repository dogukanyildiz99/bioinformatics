%Doğukan Yıldız 20360859095
data=getgenbank('L43967','SequenceOnly',true);
[H1 F1]=dimercount(data);
IPD=sum(F1,1); %bu, bizim pi sembolüyle gösterdiğimiz vektör
ETM=zeros(4,4); % estimated transition matrix
for i=1:size(ETM,1)
    for j=1:size(ETM,2)
        ETM(i,j)=F1(i,j)/IPD(i);
    end
end

%50000 nükleotitlik bir örnek sekans oluşturalım.
nukleotitler='ACGT';
rasgele=randsample(4,1,true,IPD);
%randsample fonksiyonu IPD'deki olasılıklara göre 1-4 arası bir sayı seçer
%IPD= [0.3457 0.1578 0.1591 0.3374]
%Buna göre yüzde 34 ihtimalle 1, ... yüzde 33 ihtimalle 4 seçilir
%Böylece simüle ettiğimiz sekansın ilk nükleotiti seçilmiş olur.

sekans=[rasgele];

for i=2:50000
    rasgele=randsample(4,1,true,ETM(sekans(i-1),:));
    %kalan nükleotitler hep markov zincirine göre seçilir.
    sekans=[sekans rasgele];
end

sekans2=[];
for i=1:length(sekans)
    sekans2=[sekans2 nukleotitler(sekans(i))];
end

disp(sekans2)

[H2 F2]=dimercount(sekans2);

%F1 ile F2'nin ne kadar benzediğine bakabilirsiniz.