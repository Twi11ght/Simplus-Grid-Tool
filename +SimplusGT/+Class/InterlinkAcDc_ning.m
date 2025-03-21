% Model of an ac-dc converter for inter-linking ac and dc grids.

% Author(s): Yitong Li

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
        v_od_r;
        v_oq_r;
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
                State = {'v_d','v_q','i_d','i_q','v_d_i','v_q_i','i_d_i','i_q_i','i_gd','i_gq','v_dc','i','theta'};
            else
                error('Error: Invalid ApparatusType.');
            end
        	Input = {'v_d','v_q','v','ang_r'};
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
            % wL_ac = obj.Para(2);
            W0 = obj.Para(9);
            % L_ac  = wL_ac/W0;
            % R_ac  = obj.Para(3);
            % R_dc  = obj.Para(5);

            xwLf    = 0.05;
            Rf      = 0.01;
            xwCf    = 0.02;
            xwLc    = 0.05;
            Rc      = 0.01;
            Xov     = 0;
            Rov     = 0;

            Lf = xwLf/W0;
            Cf = xwCf/W0;
            Lc = xwLc/W0;



            % Calculate paramters
            v_gd=Vg_ac;
            v_gq=0;
            i_gd=P_ac/Vg_ac;
            i_gq=-Q_ac/Vg_ac;

            v_gdq = v_gd + 1i*v_gq;
            i_gdq = i_gd + 1i*i_gq;
            v_dq = v_gdq - i_gdq*(Rc + 1i*w*Lc);
            i_cdq = v_dq*(1i*w*Cf);
            i_dq = i_gdq - i_cdq;
            e_dq  = v_dq - i_dq*(Rf + 1i*w*Lf);


            i_d = real(i_dq);
            i_q = imag(i_dq);     % Because of conjugate "i"
            i_d_i = -real(e_dq);
            i_q_i = -imag(e_dq);
            v_d_i = -i_d;
            v_q_i = -i_q;
            v_d   = real(v_dq);
            v_q   = imag(v_dq);
            i_gd  = real(i_gdq);
            i_gq  = imag(i_gdq);
            theta = xi;
            
            v = Vg_dc;      
            v_dc = Vg_dc;
            i    = P_dc/Vg_dc;
            ang_r = 0;
            
            % ??? Temp
            v_odq_r = v_dq + (Rov + 1i*Xov)*i_gdq*(-1);
            v_od_r = real(v_odq_r);
            v_oq_r = imag(v_odq_r);
            obj.v_od_r = v_od_r;
            obj.v_oq_r = v_oq_r;

            % Get equilibrium
            % State = {'v_d','v_q','i_d','i_q','v_d_i','v_q_i','i_d_i','i_q_i','i_gd','i_gq','v_dc','theta'};
        	x_e = [v_d; v_q; i_d; i_q;v_d_i;v_q_i;i_d_i;i_q_i;i_gd;i_gq;v_dc;i; theta];
        	u_e = [v_d; v_q; v; ang_r];
        end

        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
           % Get power flow
            % P_ac    = obj.PowerFlow(1);
            % Vg_ac   = obj.PowerFlow(3);
            % P_dc    = obj.PowerFlow(6);
            % Vg_dc   = obj.PowerFlow(8);
            
           	% Get parameters
            % xC_dc = obj.Para(1);
            % xwL_ac= obj.Para(2);
            % xR_ac= obj.Para(3);
            % xwL_dc= obj.Para(4);
            % xR_dc= obj.Para(5);
            % xfidq= obj.Para(6);
            % xfvdc= obj.Para(7);
            % xfpll= obj.Para(8);
            W0= obj.Para(9);
            w0= W0;
            v_dc_o = 2;
            
 
            xwLf   = 0.05;
            Rf_ac  = 0.01;
            xwCf   = 0.02;
            xwLc   = 0.05;
            Rc_dc  = 0.02;
            N1     = 1;
            R1     = 0.05;
            xfvdq  = 250;
            xfidq  = 1000;

            Lf_ac  = xwLf/W0;
            Lf_dc  = xwLf/W0;
            Rf_dc  = 0.01;
            Cf_ac  = xwCf/W0;
            Lc_ac  = xwLc/W0;
            C_dc   = 1;

            wv_ac  = xfvdq*2*pi;
            kpv_ac = Cf_ac*wv_ac;
            kiv_ac = Cf_ac*wv_ac^2/4* 20;

            wi_ac  = xfidq*2*pi;
            kpi_ac = Lf_ac*wi_ac;
            kii_ac = Lf_ac*(wi_ac^2)/4;
            
            
            % Get states
          	v_d   	= x(1);
         	v_q   	= x(2);
          	i_d  	= x(3);
            i_q 	= x(4);
            v_d_i   = x(5);
            v_q_i   = x(6);
            i_d_i   = x(7);
            i_q_i   = x(8);
            i_gd 	= x(9);
            i_gq    = x(10);
            v_dc    = x(11);
            i       = x(12);
            theta   = x(13);

            % Get input
        	v_gd    = u(1);
            v_gq    = u(2);
            v       = u(3);   %v_dc is state. How to deal this v;
            % ang_r  = u(4);
            
            % State space equations
            % dx/dt = f(x,u)
            % y     = g(x,u)
            Controll_type =0;
            if CallFlag == 1    
            % ### Call state equation: dx/dt = f(x,u)
                v_dr = obj.v_od_r;
                v_qr = obj.v_oq_r;
                if Controll_type ==1
                    error_vd = v_dr - v_d; %- (igd*Rov-igq*Xov) * (-1);
          	        error_vq = v_qr - v_q; %- (igq*Rov+igd*Xov) * (-1);
                    i_d_r = -(error_vd*kpv_ac + v_d_i);
                    i_q_r = -(error_vq*kpv_ac + v_q_i);
                    dv_d_i = error_vd*kiv_ac;
                    dv_q_i = error_vq*kiv_ac;
    
                    error_id = i_d_r-i_d;
                    error_iq = i_q_r-i_q;
                    e_d = -(error_id*kpi_ac + i_d_i);
                    e_q = -(error_iq*kpi_ac + i_q_i);
                    di_d_i = error_id*kii_ac;            
                    di_q_i = error_iq*kii_ac;
                else
                    e_d = v_dr;
                    e_q = v_qr;
                    dv_d_i=0;
                    dv_q_i=0;
                    di_d_i=0;
                    di_q_i=0;

                end

                Pr = v_dc * i;    % need -1 ?

                % Dc Lf
                di = (v - v_dc - Rf_dc*i)/Lf_dc;

                p  = (e_d*i_d + e_q*i_q)*(-1);

                w = ((Pr - p)/v_dc*R1 + (v_dc_o - 2))*N1*W0+w0; 
                dtheta = w;

                dv_dc  = (i-p/v_dc)/C_dc;

                % Lf equation
                % e_d - v_od = -(di_ld/dt*Lf + Rf*i_ld - w*Lf*i_lq)
                % e_q - v_oq = -(di_lq/dt*Lf + Rf*i_lq + w*Lf*i_ld)
                di_d = (v_gd - e_d - Rf_ac*i_d + w*Lf_ac*i_q)/Lf_ac;
                di_q = (v_gq - e_q - Rf_ac*i_q - w*Lf_ac*i_d)/Lf_ac;

                % Cf equation
                % -(i_ld - i_od) = Cf*dv_cd/dt - w*Cf*v_cq
                % -(i_lq - i_oq) = Cf*dv_cq/dt + w*Cf*v_cd
                % dv_d = (-(i_d - i_gd) + w*Cf_ac*v_q)/Cf_ac;
                % dv_q = (-(i_q - i_gq) - w*Cf_ac*v_d)/Cf_ac;
                dv_d = 0;
                dv_q = 0;

                % Lc equation
                % v_od - v_d = -(Lc*di_od/dt + Rc*i_od - w*Lc*i_oq)
                % v_oq - v_q = -(Lc*di_oq/dt + Rc*i_oq + w*Lc*i_od)
                % di_gd = (v_gd - v_d - Rc_dc*i_gd + w*Lc_ac*i_gq)/Lc_ac;
                % di_gq = (v_gq - v_q - Rc_dc*i_gq - w*Lc_ac*i_gd)/Lc_ac;
                di_gd = 0;
                di_gq = 0;


                % x_e = [v_d; v_q; i_d; i_q;v_d_i;v_q_i;i_d_i;i_q_i;i_gd;i_gq;v_dc; theta];
            	f_xu = [dv_d; dv_q; di_d; di_q;dv_d_i;dv_q_i;di_d_i;di_q_i;di_gd;di_gq;dv_dc;di; dtheta];
                Output = f_xu;
                
            elseif CallFlag == 2
          	% ### Call output equation: y = g(x,u)
                v_dr = obj.v_od_r;
                v_qr = obj.v_oq_r;
                if Controll_type ==1
                    error_vd = v_dr - v_d; %- (igd*Rov-igq*Xov) * (-1);
          	        error_vq = v_qr - v_q; %- (igq*Rov+igd*Xov) * (-1);
                    i_d_r = -(error_vd*kpv_ac + v_d_i);
                    i_q_r = -(error_vq*kpv_ac + v_q_i);
                    error_id = i_d_r-i_d;
                    error_iq = i_q_r-i_q;
                    e_d = -(error_id*kpi_ac + i_d_i);
                    e_q = -(error_iq*kpi_ac + i_q_i);
                else
                    e_d = v_dr;
                    e_q = v_qr;
                end
                Pr = v * i;    % need -1 ?
                p  = (e_d*i_d + e_q*i_q)*(-1);
                w = ((Pr - p)/v_dc*R1 + (v_dc_o-2))*N1*W0+w0; 
                g_xu = [i_d; i_q; i; w; v_dc; theta];
                Output = g_xu;
            end
        end

    end

end     % End class definition