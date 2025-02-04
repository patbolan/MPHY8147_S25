% Function to produce a 3D animation of the Magnetization vector, to
% support visualization of Bloch Eqns
% 
% Mvec: a Nx3 vector, representing Mx, My, and Mz over N timepoints
% dt: optional delay time between each update (to slow it down)
%
function animate_Mvector(Mvec, dt, bTracer)

% Use a default for dt if not provided
if nargin<2
    dt = 0.015;
end
if nargin<3
    bTracer = false;
end

clf
hold on

% Draw the vector (from origin to Mvec) at the first timepoint, save 
% handle as hM so we can adjust its appearace
origin = [0 0 0];
hM = quiver3(origin(1), origin(2), origin(3), Mvec(1,1), Mvec(1,2), Mvec(1,3), 0);
hMxy = quiver3(origin(1), origin(2), origin(3), Mvec(1,1), Mvec(1,2), 0, 0);
set(hM, 'color', 'r');
set(hM, 'LineWidth', 2);

set(hMxy, 'color', [0.5 0 0]);
set(hMxy, 'LineWidth', 1);

% Adjust view
set(gca, 'xlim', [-1 1].*1.5)
set(gca, 'ylim', [-1 1].*1.5)
set(gca, 'CameraPosition', [12 6 6]);
axis off


% Display X,Y,Z axes
axSz = 1.3;
hX = quiver3(0, 0, 0, axSz, 0, 0, 'k');
hY = quiver3(0, 0, 0, 0, axSz, 0, 'k');
hZ = quiver3(0, 0, 0, 0, 0, axSz, 'k');

% Loop over all other points and update the vector
for idx=2:size(Mvec, 1)
    
    hM.UData = Mvec(idx,1);
    hM.VData = Mvec(idx,2);
    hM.WData = Mvec(idx,3);
    
    hMxy.UData = Mvec(idx,1);
    hMxy.VData = Mvec(idx,2);

    % Tracer dot. 
    if bTracer
        Mlast = [Mvec(idx-1,1), Mvec(idx-1,2), Mvec(idx-1,3)];
        %hL = plot3([Mlast(1) Mvec(idx,1)], [Mlast(2) Mvec(idx,2)], [Mlast(3) Mvec(idx,3)], '.');
        hL = plot3(Mvec(idx,1), Mvec(idx,2), Mvec(idx,3), '.');
        set(hL, 'Color', [1 1 1].*.01);
    end
    
    pause(dt);
    
end
hold off
fprintf('Animation complete.\n');





