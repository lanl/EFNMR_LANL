% Earth's Field NMR Simulation for 2,6-difluoropyridine
% Replicates simulated spectra in Figure 7 of https://doi.org/10.1016/j.jmr.2023.107540
% Inlcudes a weighted addition of uncoupled 1H signals
%
% Adam Altenhof; adamaltenhof@gmail.com
% Derrick Kaseman; kaseman1@llnl.gov

clear all

%14N relaxation rates
% N2 = logspace(log10(1), log10(1e6), 16); %good way to sample in general
N2 = [0, 1, 10, 100, 500, 1e3, 1e4, 1e5, 1e6]; 


for i = 1:numel(N2)
    clear inter bas sys

    sys.enable={'greedy','gpu'};  

    % Sequence parameters
    parameters.sweep=1/1E-4;
    parameters.npoints=60000;
    parameters.zerofill=2^19;
    parameters.offset=0;
    parameters.spins={'1H'};
    parameters.axis_units='Hz';
    parameters.invert_axis=0;
    parameters.detection = 'uniaxial';
    parameters.flip_angle = pi/2; 
    
    
    H=2312.45;  %Field strength
    sys.magnet=H/(spin('1H')/(2*pi));
    sys.isotopes={'1H','1H','1H','19F','19F','14N'};
    
    inter.zeeman.scalar = {6.98 8.06 6.98 -70.69 -70.69 0};
    
    Ha=1; %1H 
    Hb=2; %1H
    Hc=3; %1H
    X=4;  %19F 
    Y=5;  %19F
    Z=6; %14N
    
    
    J0(1)= 0.22;  %1H and 19F relaxation
    J0(2) = N2(i); %N relaxation

    J0(3) = 0.5; %uncoupled 1H relax
    w = 0.7;     %uncoupled 1H weight

    
    TAB = 7.92; % triplet A,B;
    TBF = 8.08; % triplet B,XY; 
    TBF2 = 1.29; % triplet A,Y / C,X;
  
    for I=Ha 
        inter.coupling.scalar{I,Hb}=TAB;
        inter.coupling.scalar{I,Hc}=0.55;  
        inter.coupling.scalar{I,X}=-2.47; 
        inter.coupling.scalar{I,Y}=TBF2;
        inter.coupling.scalar{I,Z}=0.71; 
        inter.r2_rates(I)=J0(1);
    end
    
    for I=Hb 
        inter.coupling.scalar{I,Hc}=TAB;
        inter.coupling.scalar{I,X}=TBF; 
        inter.coupling.scalar{I,Y}=TBF; 
        inter.coupling.scalar{I,Z}=0;
        inter.r2_rates(I)=J0(1);
    end
    
    for I=Hc 
        inter.coupling.scalar{I,X}=TBF2;
        inter.coupling.scalar{I,Y}=-2.47; 
        inter.coupling.scalar{I,Z}=1.03; 
        inter.r2_rates(I)=J0(1);
    end
    
    for I=X
        inter.coupling.scalar{I,Y}=-12.23;
        inter.coupling.scalar{I,Z}=37.35;
        inter.r2_rates(I)=J0(1);
    end
    
    for I=Y
        inter.coupling.scaler{I,Z}=37.35;
        inter.r2_rates(I)=J0(1);
    end
    
    inter.r2_rates(Z)=J0(2);
    
    inter.coupling.scalar{Y,Y}=0;
    inter.coupling.scalar{Z,Z}=0;
    
      
    %Relaxation
    inter.relaxation={'t1_t2'};
    inter.r1_rates=inter.r2_rates;
    inter.rlx_keep='diagonal';
    inter.equilibrium='zero';
    
    % Basis set
    bas.connectivity='scalar_couplings';
    bas.space_level=4;
    bas.level=4;
    bas.formalism='sphten-liouv'; 
    bas.approximation='IK-2';
    
    
    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    
    % Simulation
    fid=liquid(spin_system,@zerofield,parameters,'labframe');
    % Fourier transform
    spectrum=fftshift(fft(fid,parameters.zerofill));
    spectrum=real(spectrum)./max(real(spectrum));
    hz=linspace(-parameters.sweep/2,parameters.sweep/2,parameters.zerofill);
    
    

    %Uncoupled H
    clear inter bas
    sys.isotopes={'1H'};
    inter.r2_rates(1)=J0(3);

    %Relaxation
    inter.relaxation={'t1_t2'};
    inter.r1_rates=inter.r2_rates;
    inter.rlx_keep='diagonal';
    inter.equilibrium='zero';
    
    % Basis set
    bas.connectivity='scalar_couplings';
    bas.space_level=1;
    bas.level=1;
    bas.formalism='sphten-liouv'; 
    bas.approximation='IK-2';
    
    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);
    
    % Simulation
    fid=liquid(spin_system,@zerofield,parameters,'labframe');
    % Fourier transform
    spectrum2=fftshift(fft(fid,parameters.zerofill));
    spectrum2=real(spectrum2)./max(real(spectrum2));

   
    spectrum = spectrum+w*spectrum2;


    plot(hz,spectrum)
    xlim([2115, 2340])

    name = strcat('2_6_FP_N_rate', num2str(J0(2)));
    savefig(gcf,strcat(name,'.fig'));
    save(strcat(name,'.mat'), 'spectrum')

    close all
end


for i = 1:numel(N2)

    name = strcat('2_6_FP_N_rate', num2str(N2(i)), '.mat');
    
    load(name)
    
    plot(hz, real(spectrum)-i*0.3)
    xlim([2115, 2340])
    xlabel('Frequency (Hz)')
    hold on

end

hold off
