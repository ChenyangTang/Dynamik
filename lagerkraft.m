function rob = lagerkraft(rob)
%Inverse Dynamik fuer Roboter rob berechnen
% Die Ergebnisse werden wiederum in der Struktur rob. gespeichert
%
% Im einzelnen werden Impuls- und Drallaenderung aller Koerper berechnet und
% ueber die Jacobi-Matrizen in die zwangsfreien Richtungen projiziert.
% Im Ergebnis wird die "linke Seite der Bewegungsgleichung"
% M*ddot_q + h berechnet und in tau_id gespeichert
% Es werden alle noetigen Groessen hier berechnet

%1. Mit Null initialisieren


lager=zeros(6,rob.N_Q);
rob.L=zeros(2,rob.N_Q);


%2. Kinematik berechnen
rob=berechne_dk_positionen(rob);
rob=berechne_dk_geschwindigkeiten(rob);
rob=berechne_dk_beschleunigungen(rob);
rob=berechne_dk_jacobis(rob);
%%
%3. Berechnung fuer alle Koerper: Impuls- und Drallaenderung
 for i=length(rob.kl):-1:1
    if i==6
    %Absolutbeschleunigung des Schwerpunkts:
    rob.kl(i).Bi_ddot_r_s = rob.kl(i).Bi_ddot_r_i+tilde(rob.kl(i).Bi_dot_omega)*rob.kl(i).Bi_r_s+...
                            tilde(rob.kl(i).Bi_omega)*tilde(rob.kl(i).Bi_omega)*rob.kl(i).Bi_r_s;
    
    %Impulsaenderung - Schwerkraft
    F_g = rob.kl(i).A_i0*rob.kl(i).m*rob.B0_g;
    dot_p = rob.kl(i).m*rob.kl(i).Bi_ddot_r_s;
    F_l= dot_p-F_g;
    
    %Drallaenderung - Moment der Schwerkraft
    M_g = tilde(rob.kl(i).Bi_r_s)*F_g;
    dot_L = rob.kl(i).I_o*rob.kl(i).Bi_dot_omega+...
            tilde(rob.kl(i).Bi_omega)*rob.kl(i).I_o*rob.kl(i).Bi_omega;
    M_l = dot_L+rob.kl(i).m*...
          tilde(rob.kl(i).Bi_r_s)*rob.kl(i).Bi_ddot_r_i-M_g;
    
    else
    %Absolutbeschleunigung des Schwerpunkts:
    rob.kl(i).Bi_ddot_r_s = rob.kl(i).Bi_ddot_r_i+tilde(rob.kl(i).Bi_dot_omega)*rob.kl(i).Bi_r_s+...
                            tilde(rob.kl(i).Bi_omega)*tilde(rob.kl(i).Bi_omega)*rob.kl(i).Bi_r_s;
    
    %Impulsaenderung - Schwerkraft
    F_g = rob.kl(i).A_i0*rob.kl(i).m*rob.B0_g;
    dot_p = rob.kl(i).m*rob.kl(i).Bi_ddot_r_s;
    F_l= dot_p-F_g-rob.kl(i+1).A_iv'*lager(1:3,i+1);
    
    %Drallaenderung - Moment der Schwerkraft
    M_g = tilde(rob.kl(i).Bi_r_s)*F_g;
    dot_L = rob.kl(i).I_o*rob.kl(i).Bi_dot_omega+...
            tilde(rob.kl(i).Bi_omega)*rob.kl(i).I_o*rob.kl(i).Bi_omega;
    M_l = dot_L+rob.kl(i).m*...
          tilde(rob.kl(i).Bi_r_s)*rob.kl(i).Bi_ddot_r_i-M_g-rob.kl(i+1).A_iv'*lager(4:6,i+1)...
          -tilde(rob.kl(i).Bv_r_vi)*rob.kl(i+1).A_iv'*lager(1:3,i+1);
    end
    lager(:,i)=[F_l;M_l];
    rob.L(1,i)=sqrt(lager(1,i)^2+lager(2,i)^2+lager(3,i)^2);
    rob.L(2,i)=sqrt(lager(4,i)^2+lager(5,i)^2+lager(6,i)^2);
    %rob.tau_id= rob.tau_id+[rob.kl(i).Bi_Jt_o;rob.kl(i).Bi_Jr]'*[F;T];
end
end
