
#include <stdlib.h>
#include <iostream>
#include <iterator>
#include <iomanip>
#include <vector>
#include <time.h>
#include <float.h>
#include <math.h>

using namespace std;


static char Usage[] = "BinMix "
                      "\n\t(input data to stdin as "
                      "\n\tcount1 majority1"
                      "\n\tcount2 majority2"
                      "\n\t...)";
                      

//---------------------------------------------------------------

double compute_log_CNK (const unsigned long N, const unsigned long K) {
    double result = 0.;
    
    for (unsigned long k = 0; k < K; k ++) {
        result += log ((N-float(k))/(K-k));
    }
    
    return result;
}

//---------------------------------------------------------------

double * get_log_rates_weights (const double * rates, const double * weights, unsigned long rate_classes) {
    double *log_rates_weights = new double [rate_classes*3];
          
    for (unsigned long j = 0; j < rate_classes; j++) {
        log_rates_weights[3*j] = log (weights[j]);
        log_rates_weights[3*j+1] = log (rates[j]);
        log_rates_weights[3*j+2] = log (1-rates[j]);
    }  
    return log_rates_weights;
}


//---------------------------------------------------------------

double compute_logL (const double * coverage, const double * majority, unsigned long data_points,
                     const double * rates, const double * weights, unsigned long rate_classes,
                     bool include_constant = false) {
      
    double logL   = 0.,
          *buffer = new double [rate_classes],
          *log_rates_weights = get_log_rates_weights (rates, weights, rate_classes);
          
          
    for (unsigned long i = 0; i < data_points; i++) {
        double max_value = -DBL_MAX,
               sum = 0.;
               
        for (unsigned long j = 0; j < rate_classes; j++) {
            buffer[j] = log_rates_weights[3*j] + majority[i]*log_rates_weights[3*j+2] + (coverage[i] - majority[i])*log_rates_weights[3*j+1];
 
            if (buffer[j] > max_value) {
                max_value = buffer[j]; 
            }
        }
        
       
        for (unsigned long j = 0; j < rate_classes; j++) {
            sum += exp (buffer[j] - max_value);
        }
        
        logL += log (sum) + max_value + (include_constant? compute_log_CNK (coverage[i], coverage[i]-majority[i]): 0.0);
    }
    
    delete [] buffer;
    delete [] log_rates_weights;
    return logL;
}

//---------------------------------------------------------------

void    compute_pij ( double * pij,
                     const double * coverage, const double * majority, unsigned long data_points,
                     const double * rates, const double * weights, unsigned long rate_classes) {
      
    double  *buffer = new double [rate_classes],
            *log_rates_weights = get_log_rates_weights (rates, weights, rate_classes);
            
    unsigned long index = 0L;
            
          
    for (unsigned long i = 0; i < data_points; i++) {
        double max_value = -DBL_MAX,
               sum = 0.;
               
        for (unsigned long j = 0; j < rate_classes; j++) {
            buffer[j] = log_rates_weights[3*j] + majority[i]*log_rates_weights[3*j+2] + (coverage[i] - majority[i])*log_rates_weights[3*j+1];
 
            if (buffer[j] > max_value) {
                max_value = buffer[j]; 
            }
        }
        
       
        for (unsigned long j = 0; j < rate_classes; j++) {
            buffer[j] = exp (buffer[j] - max_value);
            sum += buffer[j];
        }
        
        for (unsigned long j = 0; j < rate_classes; j++, index++) {
            pij [index] = buffer[j] / sum;
         }
        
    }
    
    delete [] log_rates_weights;
    delete [] buffer;
}

//---------------------------------------------------------------

void    update_values (const double * pij,
                       const double * coverage, const double * majority, unsigned long data_points,
                       double * rates, double * weights, unsigned long rate_classes) {
      
            
    for (unsigned long j = 0; j <  rate_classes; j++) {
        double sum = 0.,   
               summ = 0.,
               sumc = 0.;
               
        for (unsigned long i = 0; i < data_points; i++) {
            double thisp = pij[i*rate_classes+j];
            sum += thisp;
            summ += thisp * majority[i];
            sumc += thisp * coverage [i];
        }
        weights [j] = sum / data_points;
        if (sumc == 0.0) {
            rates   [j] = 1;
        } else {
            rates   [j] = 1.-summ / sumc;
        }
    }
}

//---------------------------------------------------------------

void     random_starting_values (double * rates, double * weights, unsigned long rate_class_count, double * prev_rates, double * prev_weights) {
    double weight_sum = 0.;
    
    for (unsigned long comp = 0; comp < rate_class_count; comp += 1) {
        if (prev_rates && prev_weights && comp < rate_class_count - 1) {
            rates [comp]  =  prev_rates[comp];
            weights[comp] =  prev_weights[comp];

        } else {
            rates [comp]  =  rand()/(double)RAND_MAX;
            weights[comp] =  rand()/(double)RAND_MAX;
        }
        weight_sum += weights[comp];
    }
    for (unsigned comp = 0; comp < rate_class_count; comp += 1) {
        weights[comp] /= weight_sum;
    }

}

//---------------------------------------------------------------

void   report_model (double final_logL, double c_aic, double * rates, double * weights, unsigned long rate_class_count) {
  cout << "{\n\"Log-L\":" << setprecision(10) << final_logL << 
    ",\n\"c-AIC\":" << c_aic <<
    ",\n\"rates\" : [";

    for (unsigned long comp = 0;  comp < rate_class_count; comp += 1) {
        cout << (comp>0?',':' ') << rates[comp];
    }
    
    cout << "],\n\"weights\": [";
    for (unsigned long comp = 0;  comp < rate_class_count; comp += 1) {
        cout << (comp>0?',':' ') << weights[comp];
    }
    cout << "]\n}\n";   
}

//---------------------------------------------------------------

double   runEM (const double * coverage, const double * majority, unsigned long data_points, 
                      double * rates, double * weights, unsigned long rate_class_count) {
                      
    double  * pij       = new double [data_points*rate_class_count];


    if (rate_class_count == 1) {
        
        double summ = 0.,
               sumc = 0.;
               
        for (unsigned long i = 0; i < data_points; i++) {
            summ += majority[i];
            sumc += coverage [i];
        }
        rates[0] = 1.-summ/sumc;
        weights[0] = 1.;
                        
    } else {
    
        double lastLogL = compute_logL (coverage, majority, data_points, rates, weights, rate_class_count);
        long max_it = 0;
        while (max_it < 100) {
            compute_pij(pij, coverage, majority, data_points, rates, weights, rate_class_count);
            update_values(pij, coverage, majority, data_points, rates, weights, rate_class_count);
            double currentLogL = compute_logL (coverage, majority, data_points, rates, weights, rate_class_count);
            if (fabs (lastLogL - currentLogL) < 1e-8) break;
            lastLogL = currentLogL;
            max_it ++;
        }
    }
    
    delete [] pij;
    return compute_logL (coverage, majority, data_points, rates, weights, rate_class_count, true);                      
                      
}

//---------------------------------------------------------------

int main (int argc, const char * argv[])
{
    if (argc != 1) {
        cerr << "Usage is `" << Usage << "'." << endl;
        return 1;
    }
    
    unsigned long rate_class_count = 1;

    
    istream_iterator<unsigned long> end_of_input;      
    istream_iterator<unsigned long> number_eater (cin);
                     
                     
    vector<double> data;
                     
    while (number_eater != end_of_input) {
         data.insert (data.end(), *number_eater++);
    }
    
    if (data.size() % 2 != 0 && data.size () < 2) {
        cerr << "Expected an even number of values (at least 2" << endl;
        return 1;
    }
    
    unsigned long data_points = data.size()/2;
    srand (time(NULL)); rand ();

    double best_aic = DBL_MAX;
    
    double   * coverage = new double [data_points],
             * majority  = new double [data_points];

    for (unsigned long i = 0; i < data_points; i++) {
        if (data.at (2*i) < data.at (2*i+1)) {
            cerr << "Coverage values must not be less than majority counts" << endl;
            return 1;
        }
        coverage[i] = data.at (2*i);
        majority[i] = data.at (2*i+1);
    }
    
    bool go_on = true;
    
    double *previous_rates = NULL,
           *previous_weights = NULL;
           
    cout << "{" << endl;       
    
    while (go_on) {
    
        double * rates       =  new double [rate_class_count],
               * weights     =  new double [rate_class_count],
               * try_rates   =  new double [rate_class_count],
               * try_weights =  new double [rate_class_count];
                       
        //cerr << "Fitting " << data_points << " observations to a binomial mixture with " << rate_class_count << " components" << endl;
        
        random_starting_values (rates, weights, rate_class_count,previous_rates,previous_weights);
        double final_logL = runEM(coverage, majority, data_points, rates, weights, rate_class_count);
        
        for (unsigned long repl = 1; repl < 50*rate_class_count; repl ++) {
             if (repl < 10) {
                random_starting_values (try_rates, try_weights, rate_class_count, previous_rates,previous_weights);
             } else {
                random_starting_values (try_rates, try_weights, rate_class_count, NULL, NULL);                
             }
             double try_logL = runEM(coverage, majority, data_points, try_rates, try_weights, rate_class_count);
             if (try_logL > final_logL) {
                final_logL = try_logL;
                for (unsigned long comp = 0;  comp < rate_class_count; comp += 1) {
                    rates [comp]  = try_rates[comp];
                    weights[comp] = try_weights[comp];
                }
              }
        }
        
        
        double parameter_count = 2*rate_class_count-1.,
               c_aic           = 2*(parameter_count*data_points/(data_points - parameter_count -1) - final_logL);
        
        if (c_aic < best_aic) {
            cout << (rate_class_count > 1? ',' : ' ') << endl << '"' << rate_class_count << '"' << " : ";
            report_model ( final_logL,  c_aic,  rates,  weights, rate_class_count);
            previous_rates   = rates;
            previous_weights = weights;
            best_aic = c_aic;
        } else {
            go_on = false;
        }
        
        rate_class_count ++;
                
        if (go_on) {
            delete [] rates;
            delete [] weights;
        }
        delete [] try_rates;
        delete [] try_weights;
    }
    
    cout << "}" << endl;       
    delete [] coverage;
    delete [] majority;
    //if (previous_rates)   delete [] previous_rates;
    //if (previous_weights) delete [] previous_weights;
    
    return 0;

}

