function [D] = computemoments(tt,H)
[~,N2]=size(H);
for ii = 1:N2
    P = abs(H(:,ii)).^2;
    %tt = t(:,ii);
    m0(ii) = trapz(tt,P);
    m1(ii) = trapz(tt,tt.*P);
    m2(ii) = trapz(tt,tt.^2.*P);
end
C = cov([m0(:) m1(:) m2(:)]);
D = [mean(m0); mean(m1); mean(m2); C(1,1); C(2,2); C(3,3);C(1,2);C(1,3);C(2,3)];
end