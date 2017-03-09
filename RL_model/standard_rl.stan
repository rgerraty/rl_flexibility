data {
	int NS;//number of subjects
	int MT;//maximum number of trials
	int NC;//number of choices (2)
	int NT[NS];//number of trials per subject
	real<lower=-1,upper=1> rew[NS,MT];//subject x trial reward, -1 for missed
	int choice[NS,MT];//chosen option, -1 for missed
	int choice_two[NS,MT];//1=chose red,0=chose blue, -1 for missed
}

transformed data{
  int N;
  N = sum (NT);
}

parameters {
  //hyperpriors on alpha distribution
  real<lower=1> a1;
  real<lower=1> a2;
  
  //hyperpriors on beta distribution
  real b_mean;
  real<lower=0> b_sd;
  
  //subject-level alpha and betas
  real<lower=0,upper=1> alpha[NS];
  vector[NS] beta;
	
}


transformed parameters{
  //subject x trials x choice Q value matrix
  real<lower=0, upper=1> Q[NS,MT,NC]; 
  
  //prediction error matrix
  real delta[NS,MT];
  
  //need to define because missing trials will recreate nan's otherwise
  Q=rep_array(0.0,NS,MT,NC);
  delta=rep_array(0.0,NS,MT);
  
  for (s in 1:NS) {
  	for (t in 1:NT[s]) {
  	  
  	  //set initial values of Q and delta on first trial
		  if(t == 1) {
		    for (c in 1:NC){
		      Q[s,t,c]=0.5;
		    }
		    delta[s,t]=0;
		  }
		    if (rew[s,t] >= 0){
		      //PE = reward-expected
		      delta[s,t]=rew[s,t]-Q[s,t,choice[s,t]];
		      
		      if (t<NT[s]){
		        //update value with alpha-weighted PE
		        Q[s,t+1,choice[s,t]]= Q[s,t,choice[s,t]] + alpha[s]*delta[s,t];
		        
		        //value of unchosen option is not updated
		        Q[s,t+1,abs(choice[s,t]-3)]=Q[s,t,abs(choice[s,t]-3)];
		        
		      }
		    } else {
		        //if no response, keep Q value and set delta to 0
		        if (t<NT[s]){
		          for (c in 1:NC){
		            Q[s,t+1,c]=Q[s,t ,c];
		          }
		        }
		        delta[s,t]=0;
		    }
		}
  }
}


model {
  //hyperpriors
  a1 ~ normal(0,5);
  a2 ~ normal(0,5);
  b_mean ~ normal (0,5);
  b_sd ~ cauchy (0,2.5);
  
  //distributions of subject effects
  alpha ~ beta(a1,a2);
  
  //for (s in 1:NS){
    beta~normal(b_mean,b_sd);
  //}
  
  
  //data generating process (likelihood)
	for (s in 1:NS) {
		for (t in 1:NT[s]) {
		  if (choice[s,t] > 0) {
		    choice_two[s,t] ~ bernoulli_logit(
		      beta[s]*(Q[s,t,2]-Q[s,t,1]));
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
        beta[s]*(Q[s,t,2]-Q[s,t,1]));
      }
      n = n+1;
    }
  }
}
