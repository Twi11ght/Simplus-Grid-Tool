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
                State = {'i_d','i_q','v_dc','i','theta'};
            else
                error('Error: Invalid ApparatusType.');
            end
        	Input = {'v_d','v_q','v','ang_r'};        % ?ang_r
            Output = {'i_d','i_q','i','w','v_dc','theta'};
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
            % P_dc    = obj.PowerFlow(6);
            P_dc = -P_ac;
            Vg_dc   = obj.PowerFlow(8);

            % Get parameters
            wL_ac = obj.Para(2);
            W0 = obj.Para(9);
            L_ac  = wL_ac/W0;
            R_ac  = obj.Para(3);
            R_dc  = 0.008;

            % Calculate paramters
            i_d = P_ac/Vg_ac;
            i_q = -Q_ac/Vg_ac;     % Because of conjugate "i"
           
            v_d = Vg_ac;
            v_q = 0;
            % syms i_d i_q
            % equ_id = (v_d - 1 + w*L_ac*i_q - R_ac*i_d)/L_ac;
            % equ_iq = (v_q - 0 - w*L_ac*i_d - R_ac*i_q)/L_ac;
            % equ_i  = [equ_id ,equ_iq];
            % [i_d,i_q] = solve(equ_i,[i_d,i_q]);
            % i_d = single(i_d);
            % i_q = single(i_q);

            theta = xi;
            ang_r =0;
            v = Vg_dc;
            i = P_dc/Vg_dc;

            v_dc = v - R_dc*i;

            % Get equilibrium
        	x_e = [i_d; i_q; theta; v_dc; i];
        	u_e = [v_d; v_q; v;ang_r];
        end

        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
            % Get power flow
            % P_ac    = obj.PowerFlow(1);
            % Vg_ac   = obj.PowerFlow(3);
            % Vg_dc   = obj.PowerFlow(8);

           	% Get parameters
            % xC_dc = 1;
            % xwL_ac= 0.05;
            % xR_ac= 0.01;
            % xwL_dc= 0.05;
            % xR_dc= 0.01;
            % xC_dc  = obj.Para(1);
            % xwL_ac = obj.Para(2);
            % xR_ac  = obj.Para(3);
            % xwL_dc = obj.Para(4);
            % xR_dc  = obj.Para(5);
            % xfidq  = obj.Para(6);
            % xfvdc  = obj.Para(7);
            % xfpll  = obj.Para(8);
            W0     = obj.Para(9);

            xwL_dc =0.08;
            xwL_ac =0.05;
            L_ac = xwL_ac/W0;
            R_ac = 0.01;
            L_dc = xwL_dc/W0;
            R_dc = 0.008;
            C_dc = 1;

            N = 1*W0;
            R = 0.05;            
            v_dc_ref = 2;
            w_ref = W0;

            % Get states
            i_d   = x(1);
            i_q   = x(2);
            theta = x(3);
            v_dc  = x(4);
            i     = x(5);

            % Get input
            v_d   = u(1);
            v_q   = u(2);
            v     = u(3);

            % State space equations
            % dx/dt = f(x,u)
            % y     = g(x,u)

            % Voltage control (open loop)
            e_d = 1;
            e_q = 0;

            % Matching control
            p1 = (e_d*i_d + e_q*i_q)*(-1);
            pr1 = v_dc*i;
            dv_dc = (pr1 - p1)/v_dc/C_dc;
            w = (v_dc + (pr1 - p1)/v_dc*R - v_dc_ref)*N + w_ref;
            dtheta = w;

            % Rectifier-side inductor
            di_d = (v_d - e_d + w*L_ac*i_q - R_ac*i_d)/L_ac;
            di_q = (v_q - e_q - w*L_ac*i_d - R_ac*i_q)/L_ac;

            % DC Line
            d_i = (v - v_dc -i*R_dc)/L_dc;

            if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)

                % Output state
                f_xu = [ di_d; di_q; dv_dc; d_i; dtheta];
                Output = f_xu;

            elseif CallFlag == 2
          	% ### Call output equation: y = g(x,u)
                g_xu = [i_d; i_q; i; w; v_dc; theta];
                Output = g_xu;
            end
        end

    end

end     % End class definition

