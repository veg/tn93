
        if (c1 < 4UL) { // c1 resolved and c2 is not
          if (matchMode != SKIP) {
            if (resolutionsCount[c2] > 0.) {
              if (matchMode == RESOLVE || (matchMode == SUBSET && resolveTheseAmbigs[c2]))
                if (resolutions[c2][c1]) {
                  ambig_count ++;
                  integer_counts[c1][c1] ++;
                  continue;
                }
              
              if (resolutions[c2][0])
                float_counts[c1][0] += resolutionsCount[c2];
              if (resolutions[c2][1])
                float_counts[c1][1] += resolutionsCount[c2];
              if (resolutions[c2][2])
                float_counts[c1][2] += resolutionsCount[c2];
              if (resolutions[c2][3])
                float_counts[c1][3] += resolutionsCount[c2];
            }
          }
        }
        else {
          if (matchMode != SKIP) {
            if (c2 < 4UL) { // c2 resolved an c1 is not
              if (resolutionsCount[c1] > 0.) {
                if (matchMode == RESOLVE || (matchMode == SUBSET && resolveTheseAmbigs[c1])) {
                  if (resolutions[c1][c2]) {
                    ambig_count ++;
                    integer_counts[c2][c2] ++;
                    continue;
                  }
                }
                
                if (resolutions[c1][0])
                  float_counts[0][c2] += resolutionsCount[c1];
                if (resolutions[c1][1])
                  float_counts[1][c2] += resolutionsCount[c1];
                if (resolutions[c1][2])
                  float_counts[2][c2] += resolutionsCount[c1];
                if (resolutions[c1][3])
                  float_counts[3][c2] += resolutionsCount[c1];
              }
            } else {
              // ambig and ambig
              double norm = resolutionsCount[c1] * resolutionsCount[c2];
              //cout << int(c1) << ":" << int(c2) << "/" << norm << endl;
              if (norm > 0.0) {
                if (matchMode == RESOLVE || (matchMode == SUBSET && resolveTheseAmbigs[c1] && resolveTheseAmbigs[c2])) {
                  ambig_count ++;
                  long matched_count = 0L,
                  positive_match [4] = {0,0,0,0};
                  for (long i = 0; i < 4L; i ++) {
                    if (resolutions[c1][i] && resolutions[c2][i]) {
                      matched_count ++;
                      positive_match[i] = 1;
                    }
                  }
                  
                  if (matched_count > 0L) {
                    double norm2 = 1./matched_count;
                    
                    for (long i = 0; i < 4L; i ++) {
                      if (positive_match[i]) {
                        float_counts[i][i] += norm2;
                      }
                    }
                    continue;
                  }
                }
                
                for (long i = 0; i < 4L; i ++) {
                  if (resolutions[c1][i]) {
                    for (long j = 0; j < 4L; j ++) {
                      if (resolutions [c2][j]) {
                        float_counts[i][j] += norm;
                      }
                    }
                  }
                }
              }
            }
          }
        }

