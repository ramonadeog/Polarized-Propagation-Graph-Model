function [A] = genWeightAdjacency(g,Ns,E,Nr,Nt,f,Rc,Psi_e)
%GENWEIGHTADJACENCY 
%Generate weighted adjacency matrix A(f) for a unipolarized graph 
%with Nr receiver(s) and Nt transmitter(s) at freq. f
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
%rng('default')
%%
N = Ns+Nr+Nt;
%E(1:Nt,:) = 0;                                %Transmitter is a source
%E(:,Nt+1:Nt+Nr) = 0;                          %Receiver is a sink
%E(Nt+1:Nt+Nr,1:Nt) = 1;                       %Include Direct path
%E(logical(eye(size(E))))= 0;                  %No loops
%
% Rscx = Rgrid(1)*rand(Ns,1);
% Rscy = Rgrid(2)*rand(Ns,1);
% Rscz = Rgrid(3)*rand(Ns,1);
% Rsc =[Rscx(:) Rscy(:) Rscz(:)];

dis_e = pdist2(Rc,Rc);  %Distance on all possible edges
tau_e = dis_e/3e8;      %Delay on all possible edges


Dtau = tau_e(Nt+1:Nt+Nr,1:Nt);
Ttau = tau_e(Nt+Nr+1:Nt+Nr+Ns,1:Nt);
Tcard = (find(Ttau));
Rtau = tau_e(Nt+1:Nt+Nr,Nt+Nr+1:Nt+Nr+Ns);
Rcard = (find(Rtau));
Btau = tau_e(Nt+Nr+1:Nt+Nr+Ns,Nt+Nr+1:Nt+Nr+Ns);   
Bcard = (find(Btau));
deg_B = zeros(Ns,Ns);
Bg = zeros(Ns,Ns);
BEdge = E(Nt+Nr+1:Nt+Nr+Ns,Nt+Nr+1:Nt+Nr+Ns);
for jj  = 1:Ns
    deg_B(:,jj) = ones(1,Ns)*length(find(Btau(:,jj)));
end
 %mus = mean(Btau(Bcard))*1e9;
 %g = db2pow(rho*mus/2);
Bg(Bcard) = (g^2.*BEdge(Bcard))./deg_B(Bcard);
Ag = zeros(N,N);
%Compute Edge gain matrix.
Dg = 1./(4*pi*f.*Dtau);
Ag(Nt+1:Nt+Nr,1:Nt) = Dg;
Tg = zeros(Ns,Nt);
Tg(Tcard) = (4*Nt*length(Tcard)./(4*pi*f*sum(sum(Ttau)))).*(1./(Ttau(Tcard).^(2).*sum(sum(Ttau(Tcard).^(-2)))));
Ag(Nt+Nr+1:Nt+Nr+Ns,1:Nt) = sqrt(Tg);
Rg = zeros(Nr,Ns);
Rg(Rcard) = (4*Nr*length(Rcard)./(4*pi*f*sum(sum(Rtau)))).*(1./(Rtau(Rcard).^(2).*sum(sum(Rtau(Rcard).^(-2)))));
Ag(Nt+1:Nt+Nr,Nt+Nr+1:Nt+Nr+Ns) = sqrt(Rg);
Ag(Nt+Nr+1:Nt+Nr+Ns,Nt+Nr+1:Nt+Nr+Ns) = sqrt(Bg);
A = E.*Ag.*exp(1j*(Psi_e-2*pi*tau_e*f));

end

