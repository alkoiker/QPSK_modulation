%QPSK MODULAZIO PROGRAMA, IKER ALKORTA ETA IGOR BOGAZ
clc;
clear all;
close all;

atala=0;
while atala~=3
    atala=menu('   AUKERATU  ', ' BIT-FLUXU BATEN QPSK MODULAZIOA EGIN ', ' QPSK SEINALEA IRAKURRI FITXATEGITIK',' IRTEN');
    switch atala
        case 1 %QPSK MODULAZI0A

            %%%%Datuen sarrera%%%
            datuak = sartuSarrerakoBitak;
         
            %Banatu bakoiti(inphase) eta bikoiti(quadrature) posizioak (MULTIPLEXAZIOA)
            [odd,even]= multiplexazioa(datuak);

            %DATUAK
            sinbolo_rate = input("Sartu sinbolo maiztasuna (sinbolo/s): ");%Transmisioaren bit rate-a
            %fc=sinbolo_rate*2; % Garraiatzailearen frekuentzia, periodoko sinbolo bat
            fc=input("Sartu garraiatzailearen frekuentzia (Hz): ");
            T=1/sinbolo_rate; % Bit bakoitzaren denbora
            lagin_bit = input("Sartu laginketa maiztasuna (Hz): ");
            t=T/lagin_bit:T/lagin_bit:T; % Bit batentzako zenbat lagin denboran
            
            %QPSK MODULAZIOA
            [QPSK_irteera, inphase_seinalea, quadrature_seinalea] = QPSKmodulazioa(datuak, even, odd, fc, t);
            
   
            %IRUDIKATU
            irudikatu(QPSK_irteera,inphase_seinalea,quadrature_seinalea, datuak, lagin_bit, T);
            
            %Seinalearen laginak CSVan idatzi
            fitxategianIdatzi(QPSK_irteera);
           

        case 2 %SEINALEA FITXATEGITIK IRAKURRI
                fitxategiaIrakurri();
                
                continue
            
        case 3 %IRTEN
            disp("Eskerrik asko")
            continue
    end
end

%*******************************************************************************************%
%BIT FLUXUA SORTZEKO FUNTZIIOA
function datuak = sartuSarrerakoBitak()
    aukera=input("Ongi etorri. Aukeratu nahi duzun bit-fluxua.\n1. Bi biteko aukera guztiak (00 01 10 11)." + ...
        "\n2. Ausazko bit-fluxua\nSartu aukera (1 edo 2): ");
    while aukera~=1 && aukera~=2
        aukera=input("Aukera okerra. 1 edo 2 zenbakiak sartu, mesedez (1 edo 2): ");
    end
    if aukera==1
        datuak=[0 0 0 1 1 0 1 1];
    elseif aukera==2
        bit_kopurua = 1;
        while(mod(bit_kopurua,2)~=0)
            bit_kopurua = input("Sartu modulatu nahi dituzun bit kopurua (bikoitia izan behar da): ");
            if mod(bit_kopurua,2)~=0
                disp('Aukera okerra. Bit kopurua bikoitia izan behar da.')
            end
        end
        for j=0:bit_kopurua
            datuak=round(rand(1,bit_kopurua));
        end
    end
end

%*******************************************************************************************%
%MULTIPLEXAZIOA EGITEKO FUNTZIOA.  bakoiti(inphase) eta
%bikoiti(quadrature), POSIZIOTAN OINARRITUZ.
function [odd, even] = multiplexazioa(datuak)
    odd= [];even=[];
    for k=1:length(datuak)
                if(mod(k,2)==0)
                    even = [even datuak(k)];
                else
                    odd = [odd datuak(k)];
                end
     end
end

%*******************************************************************************************%
%NRZ kodeketa eta QPSK modulazioa egiteko funtzioa
function [ QPSK_irteera,inphase_seinalea,quadrature_seinalea] = QPSKmodulazioa(datuak, even, odd, fc, t)
    %%%NRZ kodifikatzailea%%%%%%
    NZR_bikoitiak=2*even-1; NZR_bakoitiak=2*odd-1; % NZR kodeketa 1 denean +1 volt / 0 denean -1 volt
    
    %Modulazioa. Inphase eta quadrature tentsioak seinale garraiatzailearekin
    %biderkatuko dira, bata cosinu batekin eta bestea sinuarekin. Gero, sortu diten bi seinaleak batu egingo dira QPSK modulazioaren
    %irteerako seinalea sortzeko.
    QPSK_irteera=[];inphase_seinalea=[];quadrature_seinalea=[];

    for(i=1:length(datuak)/2)
        inphase=NZR_bakoitiak(i)*cos(2*pi*fc*t); % inphase konponentea (bit bakoitiak erabiltzen dira)
        quadrature=NZR_bikoitiak(i)*sin(2*pi*fc*t);% Quadrature konponentea (bit bikoitiak erabiltzen dira)
        inphase_seinalea=[inphase_seinalea inphase]; % inphase seinalea
        quadrature_seinalea=[quadrature_seinalea quadrature]; %quadrature seinale
        QPSK_irteera=[QPSK_irteera inphase+quadrature]; % QPSK seinalea
    end
end

%*******************************************************************************************%
%GRAFIKOAK IRUDIKATZEKO FUNTZIOA
function irudikatu(QPSK_irteera,inphase_seinalea,quadrature_seinalea, datuak, lagin_bit, T)
    tt=T/lagin_bit:T/lagin_bit:(T*length(datuak))/2;
     figure(1)
     subplot(4,1,1);
     stairs([datuak,datuak(end)],'LineWidth',2); title('Bit fluxua');

     subplot(4,1,2);
     plot(tt,inphase_seinalea,'linewidth',3), grid on;
     title('Inphase seinalea');xlabel('time(s)'); ylabel('tentsioa(V)');

     subplot(4,1,3);
     plot(tt,quadrature_seinalea,'linewidth',3), grid on;
     title(' Quadrature seinalea');xlabel('time(s)');ylabel('tentsioa(V)');

     subplot(4,1,4);
     plot(tt,QPSK_irteera,'r','linewidth',3), grid on;
     title('QPSK irteerako seinalea');xlabel('time(s)');ylabel('tentsioa(V)');
           
end

%*******************************************************************************************%
%CSV FITXATEGIAN IDAZTEKO FUNTZIOA
function fitxategianIdatzi(QPSK_irteera)
    bai ='y';ez='n';aukera=input("Seinalea .CSV fitxategi batean gorde nahi duzu? (y/n): ","s");%Fitxategia gordetzeko
    s1=strcmp(bai,aukera);
    s2=strcmp(ez,aukera);
    
    while(s1==false && s2==false)
        aukera=input("Sarrera desegokia (y/n). Seinalea .CSV fitxategi batean gorde nahi duzu? (y/n): ","s");
        s1=strcmp(bai,aukera);
        s2=strcmp(ez,aukera);
    end
        if s1==true
            q = quantizer('Float',[16 2]); QPSK_kuantifikatua = num2bin(q, QPSK_irteera);%Sortu diren QPSK laginak kuantifikatu 16 bit/lagin -etan
            fitxategia=input('Nola deitu nahi diozu fitxategiari? Sartu izena, mesedez: ','s');
            fitxategi_izena=strcat(fitxategia,'.csv');
            fileID = fopen(fitxategi_izena,'w');fprintf(fileID,'%s',QPSK_kuantifikatua);fclose(fileID);%Fitxategia idatzi.
            disp("Fitxategia gordeta.")
            
        else
            disp("Fitxategia ez da gorde.")
         end
end

%*******************************************************************************************%
%SEINALE MODULATUA IRAKURTZEKO FUNTZIOA
function fitxategiaIrakurri()
    fitxategia=input('Zer seinale ikusi nahi duzu? Sartu fitxategiaren izena, mesedez(ez jarri .csv):  ','s');
    fitxategi_izena=strcat(fitxategia,'.csv');
    if isfile(fitxategi_izena)
       % QPSK_irakurria=readmatrix(fitxategi_izena)
        fileID = fopen(fitxategi_izena,'r');formatSpec = '%s';A = fscanf(fileID,formatSpec); %Fitxategia irakurri
        QPSK_irakurria=reshape(A,[],16); %Bitak 16ka banatu

        q = quantizer('Float',[16 2]);
        QPSK_irakurria1 = bin2num(q,QPSK_irakurria);%bitak float-ara pasatu
        figure(2)
        plot(QPSK_irakurria1,'r','linewidth',3), grid on;
        title('QPSK irteerako seinalea CSVtik irakurria'); xlabel('samples'); ylabel('tentsioa(V)');
    else
        disp('Fitxategi izen hori ez da existitzen.')
    end
end