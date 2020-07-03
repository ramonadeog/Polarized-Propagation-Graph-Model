function [Rsc,Esc] = scattererPlacement(Ns,Rgrid,option)
%Place scatterer uniformly within room volume and/or walls
%   Rsc: Ns by 3 matrix of scatterer coordinates
%   Ns: Number of scatterers
%   Rgrid: room dimension vector - [L W H]
%   option: 1 for room volume placement
%   2: for wall placement with a minimum of one on each wall
%   3: wall and volume placement
%rng('default')
Rscx = Rgrid(1)*rand(Ns,1);
Rscy = Rgrid(2)*rand(Ns,1);
Rscz = Rgrid(3)*rand(Ns,1);
Esc = ones(Ns,Ns);
if option == 1
    Rsc =[Rscx(:) Rscy(:) Rscz(:)];
elseif option == 2
    if Ns < 6
        error('Number of scatter must be equal or greater than 6 for wall placement')
    end
%%%Place scatterers on walls
    indw = [randperm(6) randi([1 6],1,Ns-6)];
    indx1 = find(indw == 1);  indx2 = find(indw == 2);
    indy1 = find(indw == 3);  indy2 = find(indw == 4); 
    indz1 = find(indw == 5);  indz2 = find(indw == 6); 
    Rscx(indx1) = 0;  Rscx(indx2) = 3;
    Rscy(indy1) = 0;  Rscy(indy2) = 4;
    Rscz(indz1) = 0;  Rscz(indz2) = 3;
    Esc(indx1,indx1) = 0;
    Esc(indx2,indx2) = 0;
    Esc(indy1,indy1) = 0;
    Esc(indy2,indy2) = 0;
    Esc(indz1,indz1) = 0;
    Esc(indz2,indz2) = 0;
    Rsc =[Rscx(:) Rscy(:) Rscz(:)];
elseif option == 3
    if Ns < 7
        error('Number of scatter must be equal or greater than 7 for wall and room placement')
    end
    exW = randi([1 Ns-6]);  
    indw = [randperm(6) randi([1 6],1,exW)];
    indx1 = find(indw == 1);  indx2 = find(indw == 2);
    indy1 = find(indw == 3);  indy2 = find(indw == 4); 
    indz1 = find(indw == 5);  indz2 = find(indw == 6); 
    Rscx(indx1) = 0;  Rscx(indx2) = 3;
    Rscy(indy1) = 0;  Rscy(indy2) = 4;
    Rscz(indz1) = 0;  Rscz(indz2) = 3;
    Rsc =[Rscx(:) Rscy(:) Rscz(:)];
end
%==============Wall palcement end==============================

end

