% Computes summary statistics from graph model given parameter values

% g_min = 0.1; g_max = 0.9;
% N_min = 5; N_max = 50;
% P_min = 0.1;  P_max = 0.9;
% gamma_min = 0.0; gamma_max = 1;   %prior
% K = 2000;                 %number of parameter realization
% % modParam(1,:) = g_min+(g_max-g_min)*rand(1,K);
% % modParam(2,:) = round(N_min+(N_max-N_min)*rand(1,K));  %integer
% % modParam(3,:) = P_min+(P_max-P_min)*rand(1,K);
% % modParam(4,:) = gamma_min+(gamma_max-gamma_min)*rand(1,K);
% % 
% % modParam = [0.7 20 0.6 0.5];
% 
% modParam(1,:) = 0.7;
% modParam(2,:) = 20;
% modParam(3,:) = 0.6;
% modParam(4,:) = 0.5;
% 
% sumStats = graphModelStats(modParam, 50, 1);

function [HH, hh, specRad] = graphModelStatistics(modParam, numR, K)
    %numR = 625;    %number of realizations per parameter set
    G.roomSize = [3 4 3];
    G.freq = [58 62]*1e9;
    deltaF = 5e6;
    G.Nr = 25;
    G.Nt = 25;
    %===========================
    % numR=50;
    Rtlc = [0.90; 0.85; 2.5]; %Transmit element at origin
    Rrlc = [2.47; 1.745; 2.0];
    Rtx = Rtlc(1)*ones(1,25);
    Rty = Rtlc(2)+kron([0 0.005 0.010 0.015 0.020],ones(1,5));
    Rtz = Rtlc(3)+kron(ones(1,5),-1*[0 0.005 0.010 0.015 0.020]);
    Rrx = Rrlc(1)+kron([0 0.005 0.010 0.015 0.020],ones(1,5));
    Rry =Rrlc(2)+kron(ones(1,5),[0.020 0.015 0.010 0.005 0.00]);
    Rrz = Rrlc(3)*ones(1,25);
    Rt = [Rtx(:) Rty(:) Rtz(:)];
    Rr = [Rrx(:) Rry(:) Rrz(:)];
%     [rX1, rX2] = meshgrid(Rr(:,1),Rr(:,2));
%     rxLoc(1,:) = reshape(rX1,625,1);
%     rxLoc(2,:) = reshape(rX2,625,1);
%     rxLoc(3,:) = Rr(1,3);
%     [tX1, tX2] = meshgrid(Rt(:,1),Rt(:,2));
%     txLoc(1,:) = reshape(tX1,625,1);
%     txLoc(2,:) = reshape(tX2,625,1);
%     txLoc(3,:) = Rt(1,3);
    G.rxLoc = Rr;
    G.txLoc = Rt;
%     G.rxLoc = scattererPlacement(numR,G.roomSize,1)';
%     G.txLoc = scattererPlacement(numR,G.roomSize,1)';
    G.numPoint = (G.freq(2)-G.freq(1))/deltaF+1;
    G.Npol = 2;
    Deltat = 1/diff(G.freq); 
    Taxis = (0:G.numPoint-1)*Deltat;


    %Compute Statistics

    %handle = waitbar(0,'Initializing waitbar...');
    for ii = 1:K
        [HH,hh, specRad] = generatePGPolaNew(modParam(:,ii),G);
        %H(:,:,ii) = HH;
        for uu = 1:G.Npol
            for vv = 1: G.Npol
                Dat(uu,vv,:) = computemoments(Taxis',(squeeze(hh(:,uu,vv,:)).'));
            end
        end
%         Data(1,ii) = (Dat(1,1,1)./Dat(2,1,1)+Dat(2,2,1)./Dat(1,2,1))/2;
%         Data(2:10,ii) = mean(mean(Dat));
        
%         Data = transpose(Data);   % Transposing summary statistics matrix to be (K x nStatistics)
        %RawData(:,:,:,ii) = Dat;
        %perCom = ii/K;
        %waitbar(perCom,handle,sprintf('%d%% along...',perCom*100))
        %ii
    end

end