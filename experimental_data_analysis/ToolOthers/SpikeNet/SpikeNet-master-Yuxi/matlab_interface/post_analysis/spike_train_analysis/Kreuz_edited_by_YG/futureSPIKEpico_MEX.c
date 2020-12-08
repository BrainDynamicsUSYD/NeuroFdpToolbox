#include "mex.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
    
{
    #define bi_spike_diffs_future_t_out plhs[0]

    #define num_pairs_in prhs[0]
    #define run_isi_lengths_ruc_in prhs[1]
    #define num_trains_in prhs[2]
    #define foll_spikes_in prhs[3]
    #define cum_isi_in prhs[4]
    #define isi_in prhs[5]
    #define foll_spikes_indy_in prhs[6]
    #define udists_in prhs[7]

    int *num_pairs, *run_isi_lengths_ruc, *num_trains, *foll_spikes_indy, pac = 0, sac, trac1, trac2, M;
    double *bi_spike_diffs_future_t, *udists, *udists2, *cum_isi, *foll_spikes, *isi, logLeft, logRight;
    const mxArray *udistsPr, *udists2Pr;
    
    num_pairs = mxGetPr(num_pairs_in);
    run_isi_lengths_ruc = mxGetPr(run_isi_lengths_ruc_in);
    num_trains = mxGetPr(num_trains_in);
    foll_spikes = mxGetPr(foll_spikes_in);
    cum_isi = mxGetPr(cum_isi_in);
    isi = mxGetPr(isi_in);
    foll_spikes_indy = mxGetPr(foll_spikes_indy_in);
    
    bi_spike_diffs_future_t_out = mxCreateDoubleMatrix(0, 0, mxREAL);
    mxSetM(bi_spike_diffs_future_t_out, *num_pairs);
    mxSetN(bi_spike_diffs_future_t_out, *run_isi_lengths_ruc);
    mxSetData(bi_spike_diffs_future_t_out, mxMalloc(sizeof(double) * *num_pairs * *run_isi_lengths_ruc));
    bi_spike_diffs_future_t = mxGetPr(bi_spike_diffs_future_t_out);
    
    M  = mxGetM(udists_in);
    
    for(trac1 = 0; trac1 < *num_trains-1; ++trac1)
       for(trac2 = trac1 + 1; trac2 < *num_trains;  ++trac2) {
            pac++;
            
            udistsPr = mxGetCell(udists_in, trac2 * M + trac1);
            udists2Pr = mxGetCell(udists_in, trac1 * M + trac2);
            
            udists = mxGetPr(udistsPr);
            udists2 = mxGetPr(udists2Pr);
            
            for(sac = 0; sac < *run_isi_lengths_ruc; ++sac) {
                if (foll_spikes[trac1 + *num_trains * sac ] > foll_spikes[trac2 + *num_trains * sac]) {
                    logLeft = (foll_spikes[trac2 + *num_trains * sac] - foll_spikes[trac1 + *num_trains * sac])/2 - cum_isi[sac + 1] + foll_spikes[trac1 + *num_trains * sac];
                    logRight = (foll_spikes[trac2 + *num_trains * sac] - foll_spikes[trac1 + *num_trains * sac])/2 - cum_isi[sac] + foll_spikes[trac1 + *num_trains * sac];
                    bi_spike_diffs_future_t[(pac-1) + *num_pairs * sac] =  -(foll_spikes[trac1 + *num_trains * sac] - foll_spikes[trac2 + *num_trains * sac]
                            + udists[foll_spikes_indy[trac1 + *num_trains * sac] - 1]) / 4 *(log(logLeft + !logLeft) - log(logRight + !logRight)) / isi[sac];
                }
                else if (foll_spikes[trac1 + *num_trains * sac ] < foll_spikes[trac2 + *num_trains * sac]) {
                    logLeft = (foll_spikes[trac1 + *num_trains * sac] - foll_spikes[trac2 + *num_trains * sac])/2 - cum_isi[sac + 1] + foll_spikes[trac2 + *num_trains * sac];
                    logRight = (foll_spikes[trac1 + *num_trains * sac] - foll_spikes[trac2 + *num_trains * sac])/2 - cum_isi[sac] + foll_spikes[trac2 + *num_trains * sac];
                    bi_spike_diffs_future_t[(pac-1) + *num_pairs * sac] = -(foll_spikes[trac2 + *num_trains * sac] - foll_spikes[trac1 + *num_trains * sac]
                            + udists2[foll_spikes_indy[trac2 + *num_trains * sac] - 1]) / 4 * (log(logLeft + !logLeft) - log(logRight + !logRight)) / isi[sac];
                }
                else
                    bi_spike_diffs_future_t[(pac-1) + *num_pairs * sac] = 0;
            }
       }

    return;
}