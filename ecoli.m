%Doğukan Yıldız 20360859095
N=1000;
data=getgenbank('NC_000913','SequenceOnly',true);
sekans=data(1:N);
H=dimercount(sekans);
sayilar=cell2mat(struct2cell(H));
c=0;
E=(N-1)/16; %N-1 tane 2tuple var

for i=1:16
    if mod(i,5)==1 %AA:1, CC:6, GG:11, TT:16
        c=1+2*0.25-3*0.25^2;
    else
        c=1-3*0.25*0.25;
    end
    chisq=( sayilar(i)-E )^2 / E;
    tablo(i)=chisq/c;
end

disp(tablo')