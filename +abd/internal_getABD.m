function [ABD, ABD_INV, Qijt, NxxT, NyyT, NxyT, MxxT, MyyT, MxyT, NxxM, NyyM, NxyM, MxxM, MyyM, MxyM] =...
    internal_getABD(nPlies, Q11t, Q12t, Q16t, Q22t, Q26t, Q66t, z, nargin, deltaT, deltaM, axx, ayy, axy, bxx, byy, bxy, SECTION_POINTS)
%   Get the A, B and D matrices.
%
%   DO NOT RUN THIS FUNCTION.
%
%   Layup Analysis Tool 5.1.1 Copyright Louis Vallance 2026
%   Last modified 13-Feb-2026 11:01:37 UTC
%

%% - DO NOT EDIT BELOW LINE
%_______________________________________________________________________
%%
% Buffers for A, B and D matrix elements
Aij = zeros(3.0, 3.0); Bij = Aij; Dij = Aij;

% Buffers for thermal/moisture resultant forces and moments
NxxT = 0.0; NyyT = NxxT; NxyT = NxxT; MxxT = NxxT; MyyT = NxxT;
MxyT = NxxT; NxxM = NxxT; NyyM = NxxT; NxyM = NxxT; MxxM = NxxT;
MyyM = NxxT; MxyM = NxxT;

% Buffer for transformed reduced stiffness matrix
%{
    Qijt nominally consists of one stiffness matrix per ply (the number of
    requested sections points does not affect the value of the ABD matrix.
    However, to ensure dimensional consistency throughout the analysis,
    Qijt must contain copies of itself in each ply to account for the
    section points, which will be used later on for stress/strain
    calculation.

    For a layup with 3 section points per ply, there will be 3 copies of
    Qijt for each ply, and 3*nPlies slices of Qijt in total.
%}
Qijt = zeros(3.0, 3.0, SECTION_POINTS*nPlies);

% Index for thickness fraction buffer
iz = 2.0;

% Index for section points
spIndex = 1.0;

% Compute A, B and D matrix terms
for i = 1.0:nPlies
    % Get the themral/hydroscopic properties for the current ply
    axx_ply = axx(spIndex);
    ayy_ply = ayy(spIndex);
    axy_ply = axy(spIndex);
    bxx_ply = bxx(spIndex);
    byy_ply = byy(spIndex);
    bxy_ply = bxy(spIndex);

    % Get Qijt for the current ply 
    Qijt_ply = [Q11t(i), Q12t(i), Q16t(i);
                Q12t(i), Q22t(i), Q26t(i);
                Q16t(i), Q26t(i), Q66t(i)];

    Aij = Aij + Qijt_ply.*(z(iz) - z(iz - 1.0));
    Bij = Bij + 0.5.*Qijt_ply.*(z(iz)^2.0 - z(iz - 1.0)^2.0);
    Dij = Dij + (1.0/3.0).*Qijt_ply.*(z(iz)^3.0 - z(iz - 1.0)^3.0);

    % Compute thermal resultants
    if nargin == 5.0
        NxxT = NxxT + (deltaT*((Q11t(i)*axx_ply + Q12t(i)*ayy_ply + Q16t(i)*axy_ply)*(z(iz) - z(iz - 1.0))));
        NyyT = NyyT + (deltaT*((Q12t(i)*axx_ply + Q22t(i)*ayy_ply + Q26t(i)*axy_ply)*(z(iz) - z(iz - 1.0))));
        NxyT = NxyT + (deltaT*((Q16t(i)*axx_ply + Q26t(i)*ayy_ply + Q66t(i)*axy_ply)*(z(iz) - z(iz - 1.0))));
        MxxT = MxxT + (0.5*deltaT*((Q11t(i)*axx_ply + Q12t(i)*ayy_ply + Q16t(i)*axy_ply)*(z(iz)^2.0 - z(iz - 1.0)^2.0)));
        MyyT = MyyT + (0.5*deltaT*((Q12t(i)*axx_ply + Q22t(i)*ayy_ply + Q26t(i)*axy_ply)*(z(iz)^2.0 - z(iz - 1.0)^2.0)));
        MxyT = MxyT + (0.5*deltaT*((Q16t(i)*axx_ply + Q26t(i)*ayy_ply + Q66t(i)*axy_ply)*(z(iz)^2.0 - z(iz - 1.0)^2.0)));

        % Compute moisture resultants
        NxxM = NxxM + (deltaM*((Q11t(i)*bxx_ply + Q12t(i)*byy_ply + Q16t(i)*bxy_ply)*(z(iz) - z(iz - 1.0))));
        NyyM = NyyM + (deltaM*((Q12t(i)*bxx_ply + Q22t(i)*byy_ply + Q26t(i)*bxy_ply)*(z(iz) - z(iz - 1.0))));
        NxyM = NxyM + (deltaM*((Q16t(i)*bxx_ply + Q26t(i)*byy_ply + Q66t(i)*bxy_ply)*(z(iz) - z(iz - 1.0))));
        MxxM = MxxM + (0.5*deltaM*((Q11t(i)*bxx_ply + Q12t(i)*byy_ply + Q16t(i)*bxy_ply)*(z(iz)^2.0 - z(iz - 1.0)^2.0)));
        MyyM = MyyM + (0.5*deltaM*((Q12t(i)*bxx_ply + Q22t(i)*byy_ply + Q26t(i)*bxy_ply)*(z(iz)^2.0 - z(iz - 1.0)^2.0)));
        MxyM = MxyM + (0.5*deltaM*((Q16t(i)*bxx_ply + Q26t(i)*byy_ply + Q66t(i)*bxy_ply)*(z(iz)^2.0 - z(iz - 1.0)^2.0)));
    end

    % Update the thickness fraction index
    iz = iz + 1.0;

    % Add Qijt for the current ply over all section points for that ply
    Qijt(:, :, spIndex:spIndex + (SECTION_POINTS - 1.0)) = repmat(Qijt_ply, [1.0, 1.0, SECTION_POINTS]);

    % Update the section point loop index
    spIndex = spIndex + SECTION_POINTS;
end

% Construct the ABD matrix
ABD = [Aij, Bij; Bij, Dij];

% Get the inverse ABD matrix
ABD_INV = inv(ABD);