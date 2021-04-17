alfabe='ACGT'; %ilk 2 satır hangi baza hangi indisin
indisler=1:4; % denk geldiğini gösteriyor.

RS=['AGCT';'CTAG';'CGCG';'GGCC';'CCGG';'GATC';'TTAA';'CATG';'GTAC';'TCGA'];
%recognition site'ları tutan char matrisi
sekans=getgenbank('NC_001416','SequenceOnly',true);
N=length(sekans);
sekans_duplex=[sekans seqrcomplement(sekans)];
%kitapta olasılıklar dublex DNA üzerinden hesaplanmış
%bu nedenle seqrcomplement() fonksiyonu ile çekilen
%dizinin tümleyeni de DNA dizisinin ucuna ekleniyor.
pnuc=zeros(1,4);
%baz olasılıklarını tutacak.

sayac=1;
for i=alfabe %bu döngü baz olasılıklarını hesaplar
%bu sayılar kitapta zaten verilmiş ama biz hesaplatalım.
  pnuc(sayac)=length(find(sekans_duplex==i))/(2*N);
  sayac=sayac+1;
end

for i=1:size(RS,1) %her bir recognition site için
    H=strfind(sekans,RS(i,:)); %sekansta bu site'dan kaç tane var?
    p=1;
    for j=1:4 %her site'ın içindeki her baz için
        y=strfind(alfabe,RS(i,j)); 
        %bu site'daki j. bazın baz indisi (1-4) kaç?
        p=p*pnuc(y);
    end
    fprintf('%s\t%.5f\t%d\t%d\t%d\n',RS(i,:),p,int32(N*p),int32(N*p*(1-p)),length(H));
end