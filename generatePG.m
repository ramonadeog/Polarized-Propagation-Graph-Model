function [H,hh] = generatePG(modParam,numR,G)
%GENERATEPG 
%Generate Channel Transfer Functions from Stochastic PGM
%=======INPUTS======================================================
%modParam: 4 by K matrix
%numR: number of channel realization per model parameter set
%G: struct containing environmental/geometrical/propagation parameters 
%G.roomSize 3 by 1 vector of room dimensions
%G.freq: 2 x 1 vector of lower and upper frequency limits
%G.numPoint: number of points/samples in each transfer function
%G.rxLoc: 3 x NumR matrix of receive location
%G.txLoc: 3 x NumR matriux of transmitter location
%G.Nr: number of receive antenna location
%G.Nt: number of transmit antenna
%====OUTPUT=============================================================
%H: numR by numPoints by K
%%If you use this code or any part thereof, please consider citing our
%%paper(s)
% [1]. R. Adeogun, T. Pedersen, C. Gustafson and F. Tufvesson, "Polarimetric Wireless Indoor Channel Modeling Based on 
%Propagation Graph," in IEEE Transactions on Antennas and Propagation, vol. 67, no. 10, pp. 6585-6595, Oct. 2019.
%doi: 10.1109/TAP.2019.2925128
% [2]. R. Adeogun and T. Pedersen, "Propagation graph based model for polarized multiantenna wireless channels," 
%2018 IEEE Wireless Communications and Networking Conference (WCNC), Barcelona, 2018, pp. 1-6.
%doi: 10.1109/WCNC.2018.8377177
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Ramoni Adeogun (2019)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = linspace(G.freq(1),G.freq(2),G.numPoint);
[numPara,K] = size(modParam);
if numPara ~= 4
    modParam(4,:) = 0;
    %error('Check model parameter matrix dimension');
end
%handle = waitbar(0,'Initializing waitbar...');

for ii = 1:K
    g = modParam(1,ii); Ns = modParam(2,ii); Pvis = modParam(3,ii);
    sigma2 = modParam(4,ii);
    Nr = G.Nr; Nt = G.Nt;
    N = Ns+G.Nr+G.Nt; 
    cc=1; 
    
    while cc <= numR
        txLoc =G.txLoc(:,cc)';  
        rxLoc = G.rxLoc(:,cc)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%Adjacency Matrix Generation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        E = zeros(N,N);
        E(Nt+Nr+1:N,Nt+Nr+1:N)=(arrayfun(@(z)sum(z >= cumsum([0 1-Pvis Pvis])), rand(Ns))-1);
        E(logical(eye(size(E))))= 0;
        E(Nt+Nr+1:Nt+Nr+Ns,1:Nt) = (arrayfun(@(z)sum(z >=cumsum([0 1-Pvis Pvis])), rand(1,Ns))-1); 
        E(Nt+1:Nt+Nr,Nt+Nr+1:Nt+Nr+Ns) = (arrayfun(@(z)sum(z >=cumsum([0 1-Pvis Pvis])), rand(1,Ns))-1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Scatterer Placement Within (or on room walls) Rooms 

         Rm = scattererPlacement(Ns,G.roomSize,1); %+roomOriginLoc;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Rc = [txLoc; rxLoc; Rm];
        d = zeros(N,1);
        Psi_up = triu(0+2*pi*rand(N),1);
        Psi_e = diag(d)+Psi_up+Psi_up.';
        for ic = 1:G.numPoint
            [A] = genWeightAdjacency(g,Ns,E,Nr,Nt,f(ic),Rc,Psi_e);
            %[A] = genWeightAdjacencymu(rho,Ns,E,Nr,Nt,f(ic),Rc,Psi_e);
            B= A(Nt+Nr+1:Nt+Nr+Ns,Nt+Nr+1:Nt+Nr+Ns);
            try
              [VB] = eig(B);
            catch ME
              warning(ME.message);
              [VB] = EIG(B, 'nobalance');  
            end
            if max(abs(VB))>=1
                Hpg = 0;
                break
            end
            %%%%%%%%%%%%%%%% PG Computation=========================
            T = A(Nt+Nr+1:Nt+Nr+Ns,1:Nt); 
            R = A(Nt+1:Nt+Nr,Nt+Nr+1:Nt+Nr+Ns); %.
            Hpg(ic) = R*((eye(Ns)-B)\T); 
        end
        if Hpg == 0
            %cfail = cfail+1;
        else
            H(cc,:,ii) = Hpg;%+sigma2*sqrt(1/2)*(randn(1,G.numPoint)....
                %+1j*randn(1,G.numPoint));
            hh(cc,:,ii) = ifft(Hpg);
            cc = cc+1;      
        end   
    end
%perCom = ii/K;
%waitbar(perCom,handle,sprintf('%d%% along...',perCom*100))
end
%close(handle)
end

