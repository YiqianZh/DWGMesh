%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Digital Waveguide Mesh
% 
% Author: Chad McKell
% Date: 20 March 2020
% Place: University of California San Diego
%
% Description: This script implements a digital waveguide mesh
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
clear;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Set global parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
animate = 1;            % animate (1) or don't animate (0) the membrane
Lx = 5;                 % x-width of membrane (m)
Ly = 5;                 % y-width of membrane (m)
c = 2540;               % speed of sound in Mylar (membrane material) (m/s)
%c = 5960;               % speed of sound in Steel (membrane material) (m/s)
dur = 0.5;              % duration of simulation (seconds)
fs = 44100;             % sampling rate (Hz)
Ns = dur * fs;          % duration of simulation (samples)
nT = (0:Ns-1)/fs;       % time bin (seconds)              
Ts = 1/fs;              % sampling period (sec)
X = c * Ts;             % distance traveled in one sampling period (m)
Njx = round(Lx/X - 1);  % number of junctions along x
Njy = round(Ly/X - 1);  % number of junctions along y
xi = ceil(Njx/2);       % junction along x receiving input excitation
yi = ceil(Njy/2);       % junction along y receiving input excitation

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Define input signal
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x = [-0.5; zeros(Ns-1,1)];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Set model parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Rw1 = 0.9;              % inversion coefficient for west side of membrane
Re2 = -0.9;             % inversion coefficient for east side of membrane
Rn3 = 0.9;              % inversion coefficient for north side of membrane
Rs4 = -0.9;             % inversion coefficient for south side of membrane
w1_end = zeros(Njy,1);  % value of signal at west end of membrane
e2_end = zeros(Njy,1);  % value of signal at east end of membrane
n3_end = zeros(Njx,1);  % value of signal at north end of membrane
s4_end = zeros(Njx,1);  % value of signal at south end of membrane

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Initialize input and output ports. Note that the dimensions of each
% matrix are "number of junctions" by "number of rails".
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
portins = zeros(Njx, Njy, 4);  
portouts = zeros(Njx, Njy, 4);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute output signal
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Initialize output signal
y = zeros(1,Ns);

% Set grid for plotting
xgrid = (1:Njx)/Njx; % vector of grid points for plotting
ygrid = (1:Njy)/Njy; % vector of grid points for plotting

% Set variable to saving video
v = VideoWriter('membrane.avi');

for n = 2:Ns % sample loop
    
    for m = 1:Njx % junction loop along x
        
        for q = 1:Njy % junction loop along y
    
            % Update input ports on the UPPER HORIZONTAL rails. Note that  
            % the input ports on the upper rail are ahead one junction and   
            % ahead one time sample compared to the output ports.
            if m==1
                portins(m,q,1) = w1_end(q) * Rw1; 
            elseif m==xi && q==yi
                portins(m,q,1) = x(n-1); 
            else
                portins(m,q,1) = portouts(m-1,q,2);
            end

            % Update input ports on the LOWER HORIZONTAL rails. Note that  
            % the input ports on the lower rail are behind one junction and 
            % ahead one time sample compared to the output ports.
            if m==Njx
                portins(m,q,2) = e2_end(q) * Re2; 
            else
                portins(m,q,2) = portouts(m+1,q,1);
            end
            
            % Update input ports on the RIGHT VERTICAL rails. Note that the 
            % input ports on the right rail are ahead one junction and  
            % ahead one time sample compared to the output ports.
            if q==1
                portins(m,q,3) = n3_end(m) * Rn3; 
            else
                portins(m,q,3) = portouts(m,q-1,4);
            end

            % Update input ports on the LEFT VERTICAL rails. Note that the 
            % input ports on the left rail are behind one junction and 
            % ahead one time sample compared to the output ports.
            if q==Njy
                portins(m,q,4) = s4_end(m) * Rs4; 
            else
                portins(m,q,4) = portouts(m,q+1,3);
            end

            % Set the end-point values of the membrane. Note that these 
            % output port values correspond to one time sample in the past.
            w1_end(q) = portouts(1,q,1);
            e2_end(q) = portouts(Njx,q,2);
            n3_end(m) = portouts(m,1,3);
            s4_end(m) = portouts(m,Njy,4);
        end   
    end
    
    for m = 1:Njx % junction loop along x
        
        for q = 1:Njy % junction loop along y
   
            % Compute the pressure Pj at the 2-port junction. Note that the
            % admittance Y = S / (rho * c), where rho is the air density 
            % and c is the speed of sound in air. Since 1 / (rho * c) is 
            % constant, it can be factored out of numerator and demonator 
            % and canceled out. Furthermore, if we assume equal impedance
            % at the junctions, then S1=S2=S3=S4 and therefore S can be 
            % factored out and cancelled.
            Pj = (portins(m,q,1) + ...
                  portins(m,q,2) + ...
                  portins(m,q,3) + ...
                  portins(m,q,4))/ 2;

            % Update output ports. Note that these values will be used in  
            % the next iteration of the loop because the output ports are  
            % behind the input ports by one sample.
            portouts(m,q,1) = Pj - portins(m,q,1);
            portouts(m,q,2) = Pj - portins(m,q,2);
            portouts(m,q,3) = Pj - portins(m,q,3);
            portouts(m,q,4) = Pj - portins(m,q,4);
        end
    end
    
    if ~mod(n,1) && animate
        surf(ygrid, xgrid, portouts(:,:,1),'FaceColor','interp',...
                                           'EdgeColor','none',...
                                           'FaceLighting','gouraud','FaceAlpha',0.8);
        zlim([-0.2 0.2]);
        time = n*Ts;
        title(['t = '  num2str(time) ' sec'])
        axis off
        %view(-50,30)
        camlight left
        drawnow
        open(v);
        writeVideo(v, getframe(gca));
    end
    
    % Set output signal to the value computed at the mouth
    %y(n) = e2_end(1);
    %y(n) = portouts(ceil(Njx/4),ceil(Njy/4),1);
    y(n) = sum(sum(portouts(:,:,1)));
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Listen to output signal
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
soundsc(y,fs);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Save simulation as .wav file
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Define filename. Include parameter values in filename.
filename = 'mus206ts_mesh.wav';

% Write .wav file to MATLAB pwd at 16 bits
audiowrite(filename, y*1.5, fs, 'BitsPerSample', 16);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot time- and frequency-domain representations of y
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% subplot(2,1,1);
% plot(nT, y);
% grid on;
% title('Impulse Response');
% xlabel('Time (sec)');
% ylabel('Amplitude');
% xlim([0 0.1]);
% ylim([-1 1]);
% subplot(2,1,2);
% plotspec(y, fs,'linear');


function plotspec(x, fs, option)
% Plot the magnitude of the fft of the input x with sampling rate fs.

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute the DFT of the signal x
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N = length(x);
y = fft(x, N);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute the magnitude of the spectrum y
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y_mag = abs(y);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Normalize yMag
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y_mag = y_mag/max(y_mag);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute the frequency axis (note that this only works if the sampling 
% rate you passed in to plotspec matches the sampling rate you used to  
% define the input signal x).
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f = [0:N-1]*fs/N;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot the magnitude of the spectrum y
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if strcmp(option, 'linear') 
    plot(f, y_mag);
    xlim([0,fs/2]); % only plot to Nyquist limit
    title('Magnitude Response');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (normalized)');
    
elseif strcmp(option,'db')
    plot(f, 10*log10(y_mag));
    grid on;
    xlim([0,10000]); % fs/2]); 
    title('Filtered Signal (Frequency Domain)');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
else
    error('You must choose linear or db for the option');
end
end






