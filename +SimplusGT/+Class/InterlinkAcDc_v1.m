% Model of an FFL-VFM ac-dc converter for inter-linking ac and dc grids.

% Author(s): Yitong Li Yaoyu Hu

%% Notes
%
% The model is in 
% ac-side: load convention, admittance form.
% dc-side: load convention, admittance form.
%
% The model is simply a three phase voltage source inverter, rather than a
% back-to-back topology.

%% Class

classdef InterlinkAcDc < SimplusGT.Class.ModelAdvance
    
    % For temporary use
    properties(Access = protected)
        i_q_r;
    end
    
    methods
        % constructor
        function obj = InterlinkAcDc(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:});
        end
    end

    methods(Static)
        
        function [State,Input,Output] = SignalList(obj)
          	% Notes:
            %
            % The first three inputs must be v_d, v_q, and v; and the first
            % three outputs must be i_d, i_q, and i.
            %
            % v_d, v_q, i_d, i_q are ac-grid electrical port.v, i are the
            % dc-grid electrical port.
            %
            % v_dc is the dc link voltage, there is an inductor between v
            % and v_dc. This inductor makes the system admittance model
            % proper seen from the dc side. 
            if obj.ApparatusType==2000 || obj.ApparatusType==2001
                State = {'i_d','i_q','v_d','v_q','i_gd','i_gq','theta','v_dc','i'};
            else
                error('Error: Invalid ApparatusType.');
            end
        	Input = {'v_d','v_q','v'};        % ?ang_r
            Output = {'i_d','i_q','i','v_dc','theta','w'};
        end
        
        function [x_e,u_e,xi] = Equilibrium(obj)
            % Get the power PowerFlow values
            % The interlink converter has two sets of power flow results
            % because it is connected to two buses: the first one is ac,
            % and the second one is dc.
            P_ac    = obj.PowerFlow(1);
            Q_ac    = obj.PowerFlow(2);
            Vg_ac   = obj.PowerFlow(3);
            xi      = obj.PowerFlow(4);
            w       = obj.PowerFlow(5);
            
            P_dc    = obj.PowerFlow(6);
            Vg_dc   = obj.PowerFlow(8);
           
            % Get parameters
            wL_ac = obj.Para(2);
            W0 = obj.Para(9);
            L_ac  = wL_ac/W0;
            R_ac  = obj.Para(3);
            R_dc  = obj.Para(5);
            
            W0 = 100*pi;
            
            xwLf = 0.05;
            Rf = xwLf/5;
            Lf = xwLf/W0;

            xwCf = 0.02;
            Cf = xwCf/W0;

            xwLc = 0.5;
            Rc = xwLc/5;
            Lc = xwLc/W0;

            % Calculate paramters
            i_gd = P_ac/Vg_ac;
            i_gq = -Q_ac/Vg_ac;     % Because of conjugate "i"
            v_gd = Vg_ac;
            v_gq = 0;
            
            v_gdq = v_gd +1i*v_gq;
            i_gdq = i_gd +1i*i_gq;
            v_dq = v_gdq - i_gdq*(Rc + 1j*Lc*w);
            i_cdq = v_dq*(1i*w*Cf);
            i_ldq = i_gdq - i_cdq;

            i_d = real(i_ldq);
            i_q = imag(i_ldq);

            v_d = real(v_dq);
            v_q = imag(v_dq);

            theta = xi;
            
            v = Vg_dc;
            i = P_dc/Vg_dc;
            
            v_dc = Vg_dc - i*R_dc;
                    
            % Get equilibrium
        	x_e = [ i_d; i_q; v_d; v_q; v_dc; i_gd; i_gq; theta; i];
        	u_e = [v_gd; v_gq; v];
        end

        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get power flow
            P_ac    = obj.PowerFlow(1);
            Vg_ac   = obj.PowerFlow(3);
            Vg_dc   = obj.PowerFlow(8);
            
           	% Get parameters
            xC_dc = obj.Para(1);
            xwL_ac= obj.Para(2);
            xR_ac= obj.Para(3);
            xwL_dc= obj.Para(4);
            xR_dc= obj.Para(5);
            xfidq= obj.Para(6);
            xfvdc= obj.Para(7);
            xfpll= obj.Para(8);
            W0= obj.Para(9);
            
            W0 = 100*pi;
            xwLf = 0.05;
            xwCf = 0.02;
            xwLc = 0.5;

            Rf = xwLf/5;
            Lf = xwLf/W0;
            Cf = xwCf/W0;
            Rc = xwLc/5;
            Lc = xwLc/W0;

            ed_ref1 = 1;
            
            N = 1*W0;
            R = 0.05;
            C_dc = W0/W0;
            
            v_dc_ref = 2;
            w_ref = W0;
            
            R_dc = 0.02;
            L_dc = 0.2/W0;

            % Get states
            i_d  = x(1);
            i_q  = x(2);
            v_d  = x(3);
            v_q  = x(4);
            v_dc  = x(5);
            i_gd  = x(6);
            i_gq  = x(7);
            theta = x(8);
            i     = x(9);

            % Get input
            v_gd = u(1);
            v_gq = u(2);
            v    = u(3);

            % State space equations
            % dx/dt = f(x,u)
            % y     = g(x,u)
           
            % Voltage control (open loop)
            e_d = ed_ref1;
            e_q = 0;
                
            % Matching control
            p1 = e_d*i_d + e_q*i_q;
            pr1 = v_dc*i;
            dv_dc = (pr1 - p1)/v_dc/C_dc;
            w = (v_dc - (pr1 - p1)/v_dc*R - v_dc_ref)*N + w_ref;
            dtheta = w;

            % Rectifier-side inductor
            di_d = (e_d - v_d + w*Lf*i_q - Rf*i_d)/Lf;
            di_q = (e_q - v_q - w*Lf*i_d - Rf*i_q)/Lf;

            % Filter capacitor
            dv_d = (i_d - i_gd + w*Cf*v_q)/Cf;
            dv_q = (i_q - i_gq - w*Cf*v_d)/Cf;

            di_gd = (v_d - v_gd - i_gd*Rc + w*Lc*i_gq)/Lc;
            di_gq = (v_q - v_gq - i_gq*Rc - w*Lc*i_gd)/Lc;


            % DC Line
            d_i = (v - v_dc -i*R_dc)/L_dc;

            if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)

                % Output state
                f_xu = [ di_d; di_q; dv_d; dv_q; dv_dc; di_gd; di_gq; dtheta; d_i];
                Output = f_xu;
                
            elseif CallFlag == 2
          	% ### Call output equation: y = g(x,u)
                g_xu = [i_d; i_q; i; v_dc; theta; w];
                Output = g_xu;
            end
        end

    end

end     % End class definition