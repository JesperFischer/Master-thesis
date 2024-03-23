
from pyvirtualdisplay import Display

# Start a virtual display
with Display():
    # Import PsychoPy
    from psychopy.data import PsiHandler

import numpy as np
import math
import pandas as pd

import matplotlib.pyplot as plt


def ci(x):
    cumsum_x = np.cumsum(x) / np.sum(x)
    lower_index = np.where(cumsum_x > 0.025)[0][0]
    upper_index = np.where(cumsum_x < 0.975)[-1][-1]
    return [lower_index, upper_index]

def get_line_intervals(data, parameter):
    if parameter == "alpha":
        confidence = ci(np.mean(data, axis=1))
        rg = np.arange(-50.5, 50.5, 1)  # Assuming this is the correct equivalent for seq(0.1, 25, by = 0.1)
    elif parameter == "beta":
        confidence = ci(np.mean(data, axis=0))
        rg = np.arange(0.1, 20.1, 0.1)  # Assuming this is the correct equivalent for seq(0.1, 25, by = 0.1)
    
    
    
    upper = rg[confidence[0]]
    lower = rg[confidence[1]]
    
    data_frame = {"upper": upper, "lower": lower}
    
    return data_frame


def get_stim(lapse,alpha,beta,trials):
  df_list = []
  df = pd.DataFrame(columns=['lapse','beta','alpha','trials','X',"resp","Estimatedthreshold","Estimatedslope", "q5_threshold","q95_threshold", "q5_slope","q95_slope"])
  

  resp = np.empty(trials, dtype=object)
  Estimatedthreshold = np.empty(trials, dtype=object)
  Estimatedslope = np.empty(trials, dtype=object)
  
  Estimatedthreshold_q5 = np.empty(trials, dtype=object)
  Estimatedslope_q5 = np.empty(trials, dtype=object)

  Estimatedslope_q95 = np.empty(trials, dtype=object)
  Estimatedthreshold_q95 = np.empty(trials, dtype=object)

  f = PsiHandler(
          nTrials=trials,
          intensRange=[-50.5, 50.5],
          alphaRange=[-50.5, 50.5],
          betaRange=[0.1, 20],
          intensPrecision=1,
          alphaPrecision=1,
          betaPrecision=0.1,
          delta=0.02,
          stepType="lin",
          expectedMin=0,
      )
      
  f.next()
  for i in range(trials):
      resp[i] = np.random.binomial(1, lapse+(1-2*lapse)*(0.5+0.5*math.erf((f.intensities[i]-alpha)/(beta*np.sqrt(2)))))
      Estimatedthreshold[i] = f.estimateLambda()[0]
      Estimatedslope[i] = f.estimateLambda()[1]
      
      Estimatedthreshold_q5[i] = get_line_intervals(f._psi._probLambda[0,:,:,0],"alpha")["lower"]
      Estimatedthreshold_q95[i] = get_line_intervals(f._psi._probLambda[0,:,:,0],"alpha")["upper"]
      Estimatedslope_q5[i] = get_line_intervals(f._psi._probLambda[0,:,:,0],"beta")["lower"]
      Estimatedslope_q95[i] = get_line_intervals(f._psi._probLambda[0,:,:,0],"beta")["upper"]
      
      f.addResponse(resp[i])

      if(i != (trials-1)):
        f.next()
      else:
        break
  
    
  data = {"lapse":lapse,
          "beta":beta,
          "alpha":alpha,
          "trials":trials,
          "X":f.intensities,
          "resp":resp,
          "Estimatedthreshold":Estimatedthreshold,
          "Estimatedslope":Estimatedslope,
          "q5_threshold":Estimatedthreshold_q5,
          "q95_threshold":Estimatedthreshold_q95,
          "q5_slope":Estimatedslope_q5,
          "q95_slope":Estimatedslope_q95
          }
          

  df1 = pd.DataFrame(data)
  if not df1.empty:
    df_list.append(df1)
  


  df = pd.concat(df_list, ignore_index=True)

  

  return(df)
