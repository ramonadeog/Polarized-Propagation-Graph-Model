function [Dpol,Bpol,Tpol,Rpol] = genPolaAdjacency(Ns,Nr,Nt,Npol,G)
%GENPOLAADJENCENCY 
%A: non-weighted adjacency matrix of the graph
%Nr: Number of transmit antenna
%Npol: Number of polarization states
%G : Npol by Npol mixing matrix
%Gr: Response vectors or matrices
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
if length(G) ~= Npol
    error('size of G must be equal to Npol')
end
Gr = [1 1e-8;1e-8 1]';
Gt = [1 1e-8;1e-8 1]';
%N = Ns+Nr+Nt;
%Dpol: Polarization submtrix for direct path
Gdd = 1; %(1/sqrt(1+0.01))*[1 sqrt(0.01); sqrt(0.01) 1]; %was 0.01
Gd = Gdd.*exp(1j*(0+2*pi*zeros(Npol,Npol))); %was with rand
Dpol = kron(ones(Nr,Nt),(Gt'*Gd*Gr));   
for ii = 1:Ns
    for jj = 1:Ns
    Bpol(Npol*(ii-1)+1:ii*Npol,Npol*(jj-1)+1:jj*Npol) = G.*exp(1j*(0+2*pi*rand(Npol,Npol)));
    end
end
%Tpol: size = Ns*Npol by Nt*Npol
Tpol = zeros(Ns*Npol,Nt*Npol);
for ii = 1:Ns
    Tpol(Npol*(ii-1)+1:ii*Npol,:) = kron(ones(1,Nt),Gt'*(G.*exp(1j*(0+2*pi*rand(Npol,Npol)))));
end
Rpol = zeros(Nr*Npol,Ns*Npol);
for ii = 1:Ns
    Rpol(:,Npol*(ii-1)+1:ii*Npol) = kron(ones(Nr,1),(G.*exp(1j*(0+2*pi*zeros(Npol,Npol))))*Gr);
end
end
