%--------------------------------------------------------------------------
%
% This function initialises a density perturbation for LOKI. In this
% example, the amplitudes rho(kx,ky) are rescaled by |k| such that
% grad(E)=rho/eps0 will result in a spectrum of fields that is flat in
% k-space for EPWs.
%
% The user can specify a cutoff for the highest k mode. Modes above ~kmax/2
% should probably not be initialised.
%
% In LOKI, the user should set:
% kinetic_species.n.ic.name = "External 2D"
% kinetic_species.n.ic.file_name = "file_path" (in quotes)
% where n=species number (integer). The velocity- and temperature-dependent
% terms of the initial condition are set in the input deck.
% 
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% The USER must specify in this file the species fraction, such that
% sum_i Z_i*frac_i = 1
% i.e., multiply rho by frac_i/Z_i
% See "fract" below. If using two electron streams, that also requires
% fract=0.5 for each!
%
% In the input deck, the term kinetic_species.n.ic.frac will do nothing
% when kinetic_species.n.ic.name = "External 2D".
%
%--------------------------------------------------------------------------
%
% Annoyingly, if you want to overwrite an existing distribution, the
% current coding requires you to manually delete the old output
% distribution file.
%
%--------------------------------------------------------------------------
%
% Change log:
% - 16th October 2017: file written by Tom Chapman based on an example by
%   Bill Arrighi
%
%--------------------------------------------------------------------------

function LOKI_external_distribution_prescription

    % flag to write the hdf5 file once you're happy with rho
    do_write = 0;
    % file path of written file; should end .h5
    fout = '/p/lscratchh/chap/loki/R7016/rho_init_e2.h5';
    
    % Check rho and fft(rho) by plotting
    do_plot  = 1;

    % Mode amplitude of a k*lambda_De wave in Fourier space (note scaling
    % with k specified later)
    amp = 1.e-7;

    % nx,ny,dx,dy should agree with the input deck; do not include ghost
    % cells here
    nx = 512;
    ny = 512;

    dx = 8*pi*2/0.3/nx;
    dy = dx;
    
    kx_max = pi/dx;
    ky_max = pi/dy;
    
    % This determins where the cut-off for mode excitation in k-space is.
    % Below, this is rounded down to the nearest represented mode.
    kx_cutoff = kx_max/2.;
    ky_cutoff = ky_max/2.;
    
    % Choose carefully for ions and/or multi-species; set
    % fract = frac_i/Z_i where sum_i Z_i*frac_i = 1
    % If there are two electron species, fract=1 is probably a mistake!
    fract = 0.5;
    
    % Solution order (4 or 6); determines number of ghost cells.
    sol_order = 6;
    
%-------------------------- End of user input -----------------------------
    
    if mod(nx,2)==0
        kx =((nx-1)/nx)*linspace(-kx_max,kx_max,nx);
        dkx=mean(diff(kx));
        kx =kx-0.5*dkx;
    else
        kx =((nx-1)/nx)*linspace(-kx_max,kx_max,nx);
        dkx=mean(diff(kx));
    end
    
    if mod(ny,2)==0
        ky =((ny-1)/ny)*linspace(-ky_max,ky_max,ny);
        dky=mean(diff(ky));
        ky =ky-0.5*dky;
    else
        ky =((ny-1)/ny)*linspace(-ky_max,ky_max,ny);
        dky=mean(diff(ky));
    end
    
    x = (0:nx-1)*dx;
    y = (0:ny-1)*dy;
    x2D  = repmat(x,ny,1);
    y2D  = repmat(y.',1,nx);
    
    % centre point
    iky0 = ceil((ny+1)/2);
    ikx0 = ceil((nx+1)/2);
    
    % cut-off boundaries
    ikx_cutoffp=min(ikx0+floor(kx_cutoff/dkx),nx);
    ikx_cutoffn=max(ikx0-floor(kx_cutoff/dkx),1);
    iky_cutoffp=min(iky0+floor(ky_cutoff/dky),ny);
    iky_cutoffn=max(iky0-floor(ky_cutoff/dky),1);
    
    kx2D  = repmat(kx,ny,1);
    ky2D  = repmat(ky.',1,nx);
    
    % Set up rho in Fourier space
    rhoFFT = zeros(ny,nx);
    % Create a rectangle in Fourier space with random phases that has a
    % REAL FFT
    rhoFFT(iky0:iky_cutoffp,ikx0:ikx_cutoffp) = ...
                    exp(1i*2.*pi*rand(iky_cutoffp-iky0+1,ikx_cutoffp-ikx0+1));
    rhoFFT(iky_cutoffn:iky0,ikx_cutoffn:ikx0) = ...
                    rot90(conj(rhoFFT(iky0:iky_cutoffp,ikx0:ikx_cutoffp)),2);
    rhoFFT(iky0+1:iky_cutoffp,ikx_cutoffn:ikx0-1) = ...
                    exp(1i*2.*pi*rand(iky_cutoffp-iky0,ikx0-ikx_cutoffn));
    rhoFFT(iky_cutoffn:iky0-1,ikx0+1:ikx_cutoffp) = ...
                    rot90(conj(rhoFFT(iky0+1:iky_cutoffp,ikx_cutoffn:ikx0-1)),2);
    % Now prescribe desired mode amplitudes.
    % fftnorm gives rho an amplitude of "amp" for a single mode at given
    % kx,ky. I've paid attention to factors of 2 etc.
    % Here, rho modes are rescaled by k*lambda_De such that E fields for
    % EPWs will all have the same amplitude. This should be handled
    % differently for IAWs to reflect different scaling of E with k.
    fftnorm = 0.5*nx*ny;
    rhoFFT(iky_cutoffn:iky_cutoffp,ikx_cutoffn:ikx_cutoffp) = ...
        (amp*fftnorm) .* ...
        rhoFFT(iky_cutoffn:iky_cutoffp,ikx_cutoffn:ikx_cutoffp).* ...
        sqrt( ...
            kx2D(iky_cutoffn:iky_cutoffp,ikx_cutoffn:ikx_cutoffp).^2 ...
            +ky2D(iky_cutoffn:iky_cutoffp,ikx_cutoffn:ikx_cutoffp).^2 ...
            );
    
    % Eliminate DC term
    rhoFFT(iky0,ikx0)=0.;
    
    % Make rho; note careful use of ifft, ifftshift, and 'symmetric' - all
    % important!
    % NOTE: fract is introduced at this point
    rho = fract*(1.+ifft2(ifftshift(rhoFFT),'symmetric'));

    if do_plot
        
        figure
        imagesc([x(1),x(end)],[y(1),y(end)],rho)
        colorbar
        xlabel('x/\lambda_{De}')
        ylabel('y/\lambda_{De}')
        title('\rho')
        set(gca,'FontSize',16)
        set(gca,'YDir','normal')

        figure
        subplot(1,3,1)
        imagesc([kx(1),kx(end)],[ky(1),ky(end)],real(rhoFFT/(nx*ny)))
        colorbar
        xlabel('k_x\lambda_{De}')
        ylabel('k_y\lambda_{De}')
        title('Re[FFT(\rho)]')
        set(gca,'FontSize',16)
        set(gca,'YDir','normal')
        
        subplot(1,3,2)
        imagesc([kx(1),kx(end)],[ky(1),ky(end)],imag(rhoFFT/(nx*ny)))
        colorbar
        xlabel('k_x\lambda_{De}')
        ylabel('k_y\lambda_{De}')
        title('Im[FFT(\rho)]')
        set(gca,'FontSize',16)
        set(gca,'YDir','normal')

        subplot(1,3,3)
        imagesc([kx(1),kx(end)],[ky(1),ky(end)],abs(rhoFFT/(nx*ny)))
        colorbar
        xlabel('k_x\lambda_{De}')
        ylabel('k_y\lambda_{De}')
        title('Abs[FFT(\rho)]')
        set(gca,'FontSize',16)
        set(gca,'YDir','normal')
        
    end
    
    % note transpose of rho; Bill's convention for reading distributions
    if do_write
        
        if sol_order==4
            nghost=2;
        elseif sol_order==6
            nghost=3;
        else
            error('Unknown solution order')
        end
    
        rhoout = zeros(nx+2*nghost, ny+2*nghost);
        rhoout(nghost+1:nghost+nx,nghost+1:nghost+ny)=rho.';
        
        h5create(fout, '/2D dist', [nx+2*nghost, ny+2*nghost]);
        h5write(fout, '/2D dist', rhoout);
        fprintf('File written to: %s\n',fout);
        
    end

    fprintf('Sanity check: sum(rho)/(nx*ny) = %g\n',sum(sum(rho))/(nx*ny));
    
end
