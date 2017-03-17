data {
	int NS;//number of subjects
	int MT;//maximum number of trials
	int NStim;//number of stimuli (4)
	int NC;//number of choices (2)
	int NT[NS];//number of trials per subject
	int stim[NS,MT];//stimulus shown (1-4)
	real<lower=-1,upper=1> rew[NS,MT];//subject x trial reward, -1 for missed
	int choice[NS,MT];//chosen option, -1 for missed
	int choice_two[NS,MT];//1=chose two,0=chose one, -1 for missed
}

transformed data{
  int N;
  N = sum (NT);
}

parameters {
  //hyperpriors on alpha distribution
  real<lower=0> a1;
  real<lower=0> a2;
  
  //hyperpriors on beta distribution
  real<lower=0> b1;
  real<lower=0> b2;
  
  
  //subject-level alpha and betas
  real<lower=0,upper=1> alpha[NS];
  real<lower=0> beta[NS];
	
}


transformed parameters{
  //subject x trials x choice Q value matrix
  real Q[NS,MT,NStim,NC]; 
  
  //prediction error matrix
  real delta[NS,MT];
  
//need to assign Q and PE because missing trials will recreate nan's otherwise
  for (s in 1:NS){
    for (m in 1:MT){
      for (st in 1:NStim){
        for (c in 1:NC){
          Q[s,m,st,c]=0.0;
        }
      }
    }
  }
  
  delta=rep_array(0.0,NS,MT);
  
  for (s in 1:NS) {
  	for (t in 1:NT[s]) {
  	  
  	  //set initial values of Q and delta on first trial
		  if(t == 1) {
		    for (st in 1:NStim){
		      for (c in 1:NC){
		        Q[s,t,st,c]=.5;//want to change to 1/NC if stan ever allows int_to_real
		      }
		    }
		    delta[s,t]=0;
		  }
		    if (rew[s,t] >= 0){
		      //PE = reward-expected
		      delta[s,t]=rew[s,t]-Q[s,t,stim[s,t],choice[s,t]];
		      
		      if (t<NT[s]){
		        //update value with alpha-weighted PE
		        for (st in 1:NStim){
		          if (stim[s,t]==st){
		            Q[s,t+1,st,choice[s,t]]= Q[s,t,st,choice[s,t]] +
		            alpha[s]*delta[s,t];
		          }else{		        
		            //value of chosen option for unpresented stims not updated
		            Q[s,t+1,st,choice[s,t]]= Q[s,t,st,choice[s,t]];
		          }
		            //value of unchosen option is not updated (for any stimuli)
		            Q[s,t+1,st,abs(choice[s,t]-3)] = Q[s,t,st,abs(choice[s,t]-3)];
		      }
		    } else {
		        //if no response, keep Q value and set delta to 0
		        if (t<NT[s]){
		          for (st in 1:NStim){
		            for (c in 1:NC){
		              Q[s,t+1,st,c]=Q[s,t,st,c];
		            }
		          }
		        delta[s,t]=0;
		        }
		      }
      }
    }
  }
}


model {
  //hyperpriors
  //a1 ~ normal(0,5);
  //a2 ~ normal(0,5);
  //b_mean ~ normal (0,5);
  //b_sd ~ cauchy (0,2.5);
  
  a1 ~ cauchy(0,5);
  a2 ~ cauchy(0,5);
  b1 ~ cauchy(0,5);
  b2 ~ cauchy(0,5);
  
  //distributions of subject effects
  alpha ~ beta(a1,a2);
  beta ~ gamma(b1,1/b2);
  
  
  //data generating process (likelihood)
	for (s in 1:NS) {
		for (t in 1:NT[s]) {
		  if (choice[s,t] > 0) {
		    choice_two[s,t] ~ bernoulli_logit(
		      beta[s]*(Q[s,t,stim[s,t],2]-Q[s,t,stim[s,t],1]));
		    
		  }
		}
	}
}

generated quantities {
  real log_lik[N];
  int n;
  
  log_lik=rep_array(0,N);

  
  n = 1;
  for (s in 1:NS) {
    for (t in 1:NT[s]) {
      if (choice[s,t] > 0) {
         log_lik[n]=bernoulli_logit_lpmf(choice_two[s,t] | 
        beta[s]*(Q[s,t,stim[s,t],2]-Q[s,t,stim[s,t],1]));
      }
      n = n+1;
    }
  }
}
