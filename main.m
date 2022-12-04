clear all
close all
clc

load('DaneEksperymentu');

% Epoka 5

% Remove first 2 sec
Start5 = 123.932;
Stop5 = 151.933;
kanaly = 19;
e = 1;
fs = 500;
nazwy_kanalow={'Fp1','Fp2','F7','F3','Fz','F4','F8','T3','C3','Cz','C4','T4','T5','P3','Pz','P4','T6','O1','O2'};


% 1 wykres, 19 kanałów 

for k = 1:kanaly
    for i =1:width(dane_wynikowe.EEG_time)
        if (dane_wynikowe.EEG_time(1,i) > Start5) && (dane_wynikowe.EEG_time(1,i) < Stop5)
            Epoka_5(k,e) = dane_wynikowe.EEG_signal(k,i);
            e = e +1;
        end
    end
     e = 1;
end

    x = [1/fs:1/fs:(Stop5-Start5)];
    figure(); 
    
        for i=1:kanaly
                  subplot(4, 5, i);
                  plot(x(1,:), Epoka_5(i,:))
                  title(nazwy_kanalow(i));
                  hold on;
        end  
 %% Filtr 
 
%         rząd = 4, cz. odcięcia = 5Hz, górnoprzepustowy

  for k = 1:kanaly
[a, b] = butter(4,8/(fs/2),'high');
signalFiltered(k,:)= filter(a,b,Epoka_5(k,:));
  end

  
  
  figure(); 
  for i=1:kanaly
          subplot(4, 5, i);
          plot(x(1,:), signalFiltered(i,:))
          title("Butterworth filter - "+nazwy_kanalow(i));
          hold on;
  end 

  %% Percentyl
  
   
 for o = 1:28
  for i = 1:500
      
    EEGokna{o} = Epoka_5(:,((o-1)*500+1):1:(o*500));
  end
 end


P = [99,1,98,2,97,3,96,4,95,5,90,10];



for p = 1:2:length(P)
    
EEGokna_temp = EEGokna;


Signal_fp1 = Epoka_5(1,:);
Y = prctile(Signal_fp1,P(1,p));
    


for i=1:length(EEGokna_temp)
    if any(EEGokna_temp{i}(1,:) > Y)
        EEGokna_temp{i} = [];
        
    end
end


Signal_fp1 = Epoka_5(1,:);
Y = prctile(Signal_fp1,P(1,p+1));

for i=1:length(EEGokna_temp)
    if isempty(EEGokna_temp{i})
        continue
    elseif any(EEGokna_temp{i}(1,:) < Y)
        EEGokna_temp{i} = [];
    end
end

Signal = [];
for i=1:length(EEGokna_temp)
    Signal = [Signal EEGokna_temp{i}];
end


    % Wyświetlenie
    x_signal = [];
    x_signal = [1/fs:1/fs:length(Signal)/fs];
    figure()
        for kan=1:kanaly

                  subplot(4, 5, kan);
                  plot(x_signal(1,:), Signal(kan,:))
                  title(nazwy_kanalow(kan)+" signal <"+P(1,p+1)+";"+P(1,p)+">%");
                  hold on;
        end 

end

%% ICA

Signal_fp1 = Epoka_5;
[y, W, A] = fastica(Signal_fp1);

%%  Wyświetlenie składowych

x_signal = [1/fs:1/fs:length(y)/fs];
figure()

for kan=1:height(y)

          subplot(4, 5, kan);
          plot(x_signal(1,:), y(kan,:))
          title("ICA - składowa "+kan);
          hold on;
end 
%% Filtr
moc_y = [];
for i=1:height(y)
    
   [a, b] = butter(4,[2 4]/(fs/2),'bandpass');
   signalFiltered = filter(a,b,y(i,:));
   moc_y =[moc_y mean(signalFiltered.^2)];
   
end
    
moc_y = moc_y'

mean_moc_y = mean(moc_y);
std_y = std(moc_y);

y_po_ICA = y;

for i=1:height(y_po_ICA)
    
   if  (moc_y(i,1) > (mean_moc_y + (1.5*std_y)))
       y_po_ICA(i,:) = 0;
   end
       
end

% y -macierz składowych A-macierz demiksująca W-macierz miksująca
% moc składowej > średniej mocy ze wszystkich składowych + 1.5 std ze wszystkich składowych
%%  Wyświetlenie składowych po oczyszczeniu

x_signal = [1/fs:1/fs:length(y_po_ICA)/fs];
figure()

for kan=1:height(y_po_ICA)

          subplot(4, 5, kan);
          plot(x_signal(1,:), y_po_ICA(kan,:))
          title("ICA - składowa "+kan);
          hold on;
end 
%% Odtwarzanie sygnału

signal = W * y_po_ICA;

x = [1/fs:1/fs:length(y_po_ICA)/fs];

figure('Name','New EEG signal')
for i=1:height(signal)
          subplot(4, 5, i);
          plot(x, signal(i,:));
          nazwy_ICA = [1:height(y)];
          title(nazwy_kanalow(i));
          hold on;
end   
%% Wspólna referencja

ref = [];

for i =1: length(Epoka_5)
    ref(1,i) = mean(Epoka_5(:,i));
end

ref = ref /19;

result = [];

for i =1:length(Epoka_5)
    
    result(:,i) = (Epoka_5(:,i)) - ref(:,i);
end

    %% Wyświetlenie wspólnej referencji
    
figure()
    for kan=1:kanaly

              subplot(4, 5, kan);
              plot(x_signal(1,:), result(kan,:))
              title(nazwy_kanalow(kan)+" referencja");
              hold on;
    end 
    