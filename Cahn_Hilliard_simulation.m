% CAHN-HILLIARD SIMULATION
global u_mean u_amplitude Nx Ny Ntau tau h number_of_frames Path1 Name1 Path2 Name2 time_of_constant Tc T;
clear data1 struct1

%-----------DEFINE INITIAL PARAMETERS------------------
L=0.1;          %system length
Nx = 248;       %x lattice of the simulation
Ny = 248;       %y lattice of the simulation
Ntau_change = 20300; %number of time steps, where the system will evolve
tau = 0.09;     %time step
number_of_frames = 5;  %frame step for the display of the result
h = 1.51; %L/Nx;     %space step
u_amplitude=0.1;    %amplitude of initial concentration variation
step=1;    %step of structure averaging

%---Define paths for saving 
Path1='Yourpath\real_space';
Name1='u_2';
Path2='Yourpath\structure';
Name2='structure_2';

%---Definition of changes in the temperature. The classical Cahn-Hilliard
%assumes that T is constant all the time. In this example, there is a
%user-defined evolution.
time_of_constant=0;
Ntau=time_of_constant+Ntau_change;
Tc=1;
T_step1=0.98;
T_step2=0.05;
T_step3=0.05;
for i=1:3000
    T(i)=T_step1;
end
for i=3001:4000
    T(i)=T_step1+(T_step2-T_step1)*(i-3001)/(4000-3001);
end
for i=4001:Ntau
    T(i)=T_step3;
end
%------------------------------------------------------------



%-----------------------CALCULATION OF CAHN-HILLIARD

for i=1%6     %if different initial mean concentration values are wanted to be used, make the cycle with the corresponding change in u_mean
    u_mean=0.4;     %initial mean concentration
    disp(u_mean);
    
    % u0 = init_u();
    u_solution = zeros(Nx, Ny, Ntau);   %create an empty array of concentration
    
    for kk=1:(time_of_constant+1)  %define initial concentration distribution
        u0_simple = u_mean + (rand(Nx, Ny) - 0.5) * u_amplitude;
        %u0_sin = init_u();
        u_solution(:, :, kk) = u0_simple;
    end
    
    if Ntau<2301            %the size defines the maximal possible array to be saved in .mat on the specific computer. This threshold can be changes by user
        data1 = eulerSolver(u_solution); %solve numerically the differential equation with Euler solver. Result is the concentration map evolution with time
        struct1 = Structure(data1);   %calculate the scattering from this concentration map
        %data2 = CrankNicolsonSolver(u_solution);  %other solverscan be
        %used
    else
        step=round(Ntau/2000);  %in case the size is too big for saving as .mat, the averaging is performed
        data1=Averaging(eulerSolver(u_solution),step);
        struct1=Structure(data1);
        
    end
    
    close all;
    showData(data1, 'Euler');
    %showData(data2, 'Crank - Nicolson');
    
    
    

    %save data  
    if ~exist([Path1],'dir')
        mkdir([Path1])
    end
    save([Path1,'\',Name1,'_',num2str(u_mean),'_',num2str(u_amplitude),'_',num2str(tau),'_',num2str(Ntau),'_',num2str(T_step1),'_',num2str(T_step2),'_',num2str(T_step3),'_',num2str(h),'_wing.mat'],'data1','u_mean','u_amplitude','tau','Ntau','step','Tc','T','h')
    
    if ~exist([Path2],'dir')
        mkdir([Path2])
    end
    save([Path2,'\',Name2,'_',num2str(u_mean),'_',num2str(u_amplitude),'_',num2str(tau),'_',num2str(Ntau),'_',num2str(T_step1),'_',num2str(T_step2),'_',num2str(T_step3),'_',num2str(h),'_wing.mat'],'struct1','tau','Ntau','step','Tc','T')
    writevideo(['',Name1,'_',num2str(u_mean),'_',num2str(u_amplitude),'_',num2str(tau),'_',num2str(Ntau),'_',num2str(T_step1),'_',num2str(T_step2),'_',num2str(T_step3),'_',num2str(h),'_wing.avi'],data1)
    
end



%--------------FUNCTIONS--------------------


function showData(data, tl)
%display evolution of concentration maps
global Ntau number_of_frames Path1 Name1 tau u_mean u_amplitude;
delta = round(size(data,3)/ number_of_frames);
for z=[1:delta:size(data,3) size(data,3)]
    if ~exist([Path1,'\images'],'dir')
        mkdir([Path1,'\images'])
    end
    figure, imagesc(data(:,:,z))
    if (z==1)
        title(strcat(tl, num2str(z*tau)));
    else
        title(strcat(tl, num2str(z*tau*round(Ntau/2000))));
    end
    colormap jet
    colorbar
    set(gca,'fontsize', 20)
    set(gca,'linewidth', 2)
    pbaspect([1 1 1])
    caxis([-1 1])
    saveas(gcf,[Path1,'\images\',Name1,'_',num2str(u_mean),'_',num2str(u_amplitude),'_',num2str(tau),'_',num2str(Ntau),'_',num2str(z),'_time.png'])
    saveas(gcf,[Path1,'\images\',Name1,'_',num2str(u_mean),'_',num2str(u_amplitude),'_',num2str(tau),'_',num2str(Ntau),'_',num2str(z),'_time.fig'])
end
end

function data = CrankNicolsonSolver(u_solution)
%solve numerically differential scheme with CrankNicolsonSolver
global Ntau tau;
num_of_iterations = 3;

u_solution(:, :, 2 ) = u_solution(:, :, 1 ) + F(u_solution(:, :, 1 ))* tau; % Euler one step method

for z=3:Ntau
    %disp(z);
    u0 = u_solution(:, :, z - 1);
    u_temporary = u0;
    f0 = F(u0);
    for it = 1:num_of_iterations
        u_temporary = u0 + 1./2 * tau * ( F(u_temporary) + f0 );
    end
    u_solution(:, :, z) = u_temporary;
end

data = u_solution;

end


function data = AdamsBashforthSolver(u_solution)
%solve scheme with Adams-Bashforth Solver
global Ntau tau;

u_solution(:, :, 2 ) = u_solution(:, :, 1 ) + F(u_solution(:, :, 1 ))* tau; % Euler one step method

for z=3:Ntau
    disp(z);
    u0 = u_solution(:, :, z - 2 );
    u1 = u_solution(:, :, z - 1 );
    f1 =  F(u1);
    f0 =  F(u0);
    u_solution(:, :, z) = u1 + 1/2 * tau * (3  * f1 - f0 );
end

data = u_solution;

end

function data = eulerSolver(u_solution)
%solve scheme with euler Solver
global Ntau tau time_of_constant;
for z=(time_of_constant+2):Ntau
    disp(z);
    u0 = u_solution(:, :, z - 1 );
    u_solution(:, :, z) = u0 + F(u0,z) * tau;
end
data = u_solution;

end


function mm = m(u)   %assumption of constant mobility
%define mobility
mm = 1;
end

function kk = k(u)
%energy gradient coefficient
kk = 1;
end

function mmu =  mu(u, ddu,time)
%Cahn-Hilliard equation
global h Tc T;
mmu = -(Tc-T(time))/Tc.*u+u.^3 - 1/h^2 * ddu .* k(u);
end



function f = F(u,time)
%discretization 
global Nx Ny h;
ddu = zeros(Nx, Ny);
f = zeros(Nx, Ny);

for i=1:(Nx)
    for j=1:(Ny)
        ip1=1+mod(i+1-1,Nx);
        im1=1+mod(i-1-1,Nx);
        jp1=1+mod(j+1-1,Ny);
        jm1=1+mod(j-1-1,Ny);
        ddu(i, j) = u(ip1, j) + u(im1,j) + u(i, jm1) + u(i, jp1) - 4*u(i,j);
    end
end

mu_matrix = mu(u, ddu,time);

for i=1:(Nx)
    for j=1:(Ny)
        
        ip1=1+mod(i+1-1,Nx);
        im1=1+mod(i-1-1,Nx);
        jp1=1+mod(j+1-1,Ny);
        jm1=1+mod(j-1-1,Ny);
        
        mip = m( 1/2*(u(i,j) + u(ip1,j)) );
        mim = m( 1/2*(u(i,j) + u(im1,j)) );
        mjp = m( 1/2*(u(i,j) + u(i,jp1)) );
        mjm = m( 1/2*(u(i,j) + u(i,jm1)) );
        f(i,j) = 1/h^2 * (...
            mip * (mu_matrix(ip1, j) - mu_matrix(i, j) ) + ...
            mim * (mu_matrix(im1, j) - mu_matrix(i, j) ) + ...
            mjp * (mu_matrix(i, jp1) - mu_matrix(i, j) ) + ...
            mjm * (mu_matrix(i, jm1) - mu_matrix(i, j) )...
            );
    end
end

end


function u0 = init_u()
%periodic boundary conditions
global u_mean u_amplitude Nx Ny L V;
u0=zeros(Nx,Ny);
random=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        for kk=-V:1:V
            for l=-V:1:V
                random(i,j)=rand(1,1);
                u0(i,j) =u0(i,j)+u_amplitude*sin(2*pi/L*(kk*i*L/Nx+l*j*L/Ny+random(i,j)));
            end
        end
    end
end
u0=u0+u_mean*ones(Nx,Ny)-u_amplitude*sin(2*pi/L*(random));
end



function S = Structure(data_real)
%calculate scattering as fourier transform from the concentration map
global Nx Ny u_mean;

S = zeros(Nx-1, Ny-1, size(data_real,3));

for t=1:size(data_real,3)
    
    
    current_fft = fft2(squeeze(data_real(:,:,t)-u_mean));
    current_fft(1,:) = [];
    current_fft(:,1) = [];
    S(:,:,t)=abs(fftshift(current_fft)).^2;
    
end

end

function av=Averaging(data_av,step)
%average data with step
global Nx Ny Ntau
av=zeros(Nx,Ny,round(Ntau/step)+1);
av=data_av(:,:,1:step:Ntau);
end

function writevideo(Name,data)
%write video of the evolution of the map
global Ntau tau
figure
v=VideoWriter(Name);
open(v);
for time=1:20:size(data,3)
    imagesc(data(:,:,time))
    if (time==1)
        title(strcat('Euler', num2str(time*tau)));
    else
        title(strcat('Euler', num2str(time*tau*round(Ntau/2000))));
    end
    colormap jet
    colorbar
    set(gca,'fontsize', 20)
    set(gca,'linewidth', 2)
    pbaspect([1 1 1])
    caxis([-1 1])
    frame=getframe(gcf);
    writeVideo(v,frame);
end
end
