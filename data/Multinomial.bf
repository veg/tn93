log_factorials           = {"0": 0};

function log_fact (N) {
    if (N < Abs (log_factorials)) {
        return log_factorials[N];
    }
    
    lastN = log_factorials[Abs (log_factorials)-1];
    for (k = Abs (log_factorials); k <= N; k += 1) {
        log_factorials[k] = Log (k) + log_factorials[k-1];
    }
    
    return log_factorials[N];
    
}



function make_global_var (var_id, val) {
    ExecuteCommands ("global `var_id` = val; `var_id`:<0.999999; `var_id`:>1e-6;");
}

function precompute_mc (coverages) {
    
    sites = Rows (coverages);
    
    for (s = 0; s < sites; s+=1) {
         multinomial_coefficients[s] = log_fact(coverages[s][1]) -
                                       log_fact(coverages[s][2]) -
                                       log_fact(coverages[s][3]) -
                                       log_fact(coverages[s][4]) - 
                                       log_fact(coverages[s][5]);  
    }

    return 0;
}

function multinomial_logL (rates, weights) {
    log_rates = Transpose(Log (rates));
    logL = counts * log_rates;
    result = {sites,1};
    for (k = 0; k < sites; k+=1) {
        thisRow = logL[k][-1];
        minVal = (Min (thisRow*(-1), 1))[0];
        thisRow += minVal;
        result[k] = Log(+(thisRow["Exp(_MATRIX_ELEMENT_VALUE_)"]$weights)) - minVal + multinomial_coefficients[k];
    }
    return +result;
}

function report_site (post,rates,s,stencil) {
    result = {4,1};
    site_post = post[s][-1];
    for (k = 0; k < 4; k+=1) {
        result[k] = +(site_post $ stencil [k][-1]);
    }
    return result;
}

function posterior (rates, weights) {
    log_rates = Transpose(Log (rates));
    logL = counts * log_rates;
    rate_count = Rows (rates);
    
    result = {sites, rate_count};
    
    for (k = 0; k < sites; k+=1) {
        thisRow = logL[k][-1];
        minVal = (Min (thisRow*(-1), 1))[0];
        thisRow += minVal;
        thisRow = thisRow["Exp(_MATRIX_ELEMENT_VALUE_)"]$weights;
        thisRow = thisRow*(1/(+thisRow));
        for (i = 0; i < rate_count; i+=1) {
            result [k][i] = thisRow[i];
        }
    }
    return result;
}

LoadFunctionLibrary ("ReadDelimitedFiles");
LoadFunctionLibrary ("GrabBag");


coverage_data = ReadTabTable ("BB_full.txt", 0);
sites         = Rows         (coverage_data);

multinomial_coefficients = {sites,1};
counts = coverage_data [{{0,2}}][{{sites-1,5}}];

precompute_mc (coverage_data);

per_class = 3;
total     = per_class^4;

rates         = {total,4};
current_class = 0;

weights       = {1, total};

//generate_gdd_freqs (numberOfRates, freqs&, lfMixing&, probPrefix, incrementFlag)
Af = {}; l = "";
generate_gdd_freqs (per_class, "Af", "l", "fA", 0);
Cf = {}; 
generate_gdd_freqs (per_class, "Cf", "l", "fC", 0);
Gf = {}; 
generate_gdd_freqs (per_class, "Gf", "l", "fG", 0);
Tf = {}; 
generate_gdd_freqs (per_class, "Tf", "l", "fT", 0);


init_values = {{0.01,0.05,0.98}};

for (i = 0; i < per_class; i+=1) {
    make_global_var ("A" + i, init_values[i]);
    make_global_var ("C" + i, init_values[i]);
    make_global_var ("G" + i, init_values[i]);
    make_global_var ("T" + i, init_values[i]);
}

T0 := A0;
G0 := A0;
C0 := A0;
stencil = {4, total};

for (A = 0; A < per_class; A+=1) {
    for (C = 0; C < per_class; C+=1) {
        for (G = 0; G < per_class; G+=1) {
            for (T = 0; T < per_class; T+=1) {
                normalizer = "(A" + A + "+C" + C + "+G" + G + "+T" + T + ")";
                ExecuteCommands ("rates[current_class][0] := A" + A + "/`normalizer`");
                ExecuteCommands ("rates[current_class][1] := C" + C + "/`normalizer`");
                ExecuteCommands ("rates[current_class][2] := G" + G + "/`normalizer`");
                ExecuteCommands ("rates[current_class][3] := T" + T + "/`normalizer`");
                ExecuteCommands ("weights[current_class] := " + Af[A] + "*" + Cf[C] + "*" + Gf[G] + "*" + Tf[T]);
                stencil [0][current_class] = A > 0;
                stencil [1][current_class] = C > 0;
                stencil [2][current_class] = G > 0;
                stencil [3][current_class] = T > 0;
                current_class += 1;
            }
        }
    }
}
      
//VERBOSITY_LEVEL = 1;
Optimize (res,multinomial_logL (rates,weights));



/*for (i = 0; i < per_class; i+=1) {
    fprintf (stdout, "A", i , " = ", Eval ("A" + i), " (p = ", Eval (Af[i]), ")\n");
    fprintf (stdout, "C", i , " = ", Eval ("C" + i), " (p = ", Eval (Cf[i]), ")\n");
    fprintf (stdout, "G", i , " = ", Eval ("G" + i), " (p = ", Eval (Gf[i]), ")\n");
    fprintf (stdout, "T", i , " = ", Eval ("T" + i), " (p = ", Eval (Tf[i]), ")\n");
}*/

post = posterior (rates, weights);

fprintf (stdout, "{\"background\":", A0, ",\n\"posteriors\": [");

for (site = 0; site < sites; site +=1) {
   pp = report_site (post, rates, site, stencil);
   if (site) {
    fprintf (stdout, ",\n");
   }
   fprintf (stdout, "[", Join(",", pp), "]");
}
fprintf (stdout, "\n]}\n");