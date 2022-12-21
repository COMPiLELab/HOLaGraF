function wrapper(G::NiceGraph, h0, α_st, α_fin, thrs::Thresh, h_ε, logSizes, logSteps, logTrack, logE, logΛ; mask=ones( size(G.e) ) )
      while true
              (G.eps0 > 5.0) && break;

              if G.eps0 == h_ε
                      @printf "             alpha  stage started: \n";
                      G, thrs, track, p, Hs, Es, Λs, SZs, TrackS = alphaLevel(G, h0, α_st, α_fin, thrs; initial=true, mask=mask);
                      @printf " eps:  %f       ||||    functional :  %f  \n" G.eps0 track[end];
                      @printf "              alpha stage ended \n \n";
              else
                      G , thrs, track, p, Hs, Es, Λs, SZs, TrackS = alphaLevel(G, h0, α_st, α_fin, thrs, mask=mask); 
                      @printf " <contrained>    eps = %f    ||||  functional :  %f  || %d steps \n" G.eps0 track[end] size(track, 1);
              end
             #G.eps0 += h_ε;
             logSizes = [ logSizes; SZs ]; logSteps = [logSteps; Hs]; logTrack = [logTrack; TrackS];
              logE = [logE Es]; logΛ = [logΛ Λs];

             (track[end] < 1e-5) && break; 

              ε_target = G.eps0 + h_ε; 
              ans, track, h_log, G, L1up, L1, L0, p, e_log, λ_log  = freeGradientTransition(ε_target, G, h0, thrs, p; mask=mask);
              G.e=ans;
              @printf  "  <free_grad>  eps_fin = %f  / %f        |||| functional:   %f      || steps: %d  \n"  G.eps0 ε_target track[end] size(track, 1);
             
              logSizes = [ logSizes; size(h_log, 1) ]; logSteps = [logSteps; h_log]; logTrack = [logTrack; track];
              logE = [logE e_log]; logΛ = [logΛ λ_log];
      end
      return G, thrs, logSizes, logSteps, logTrack, logE, logΛ
end

