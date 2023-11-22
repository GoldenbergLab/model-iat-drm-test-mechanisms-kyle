functions {
    real race_pdf(real t, real alpha, real v){
        real pdf;
        pdf = alpha/sqrt(2 * pi() * pow(t, 3)) * exp(-pow(v*t-alpha, 2) / (2*t));
        return pdf;
    }
    
    real race_cdf(real t, real alpha, real v){
        real cdf;
        cdf = Phi((v*t-alpha)/sqrt(t)) + exp(2*v*alpha) * Phi(-(v*t+alpha)/sqrt(t));
        return cdf;
    }

    real race_reduce(real[] RT, int start, int end, int[] choice, int[] N_ind, int[] N_sub_ind, int[] condition, real[] mu_alpha, real[] mu_v_lower, real[] mu_v_upper, real mu_tau, real[,] alpha, real[,] v_lower, real[,] v_upper, real[] tau){

        real t;
        vector[end-start+1] prob;
        real cdf;
        real pdf;
        real out;
        int counter;
        real alpha_i;
        real v_upper_i;
        real v_lower_i;
        real tau_i;
        counter = 1;

        for (i in start:end) {
            if (N_sub_ind[N_ind[i]] != 0) {
                alpha_i = alpha[N_sub_ind[N_ind[i]], condition[i]];
                v_lower_i = v_lower[N_sub_ind[N_ind[i]], condition[i]];
                v_upper_i = v_upper[N_sub_ind[N_ind[i]], condition[i]];
                tau_i = tau[N_sub_ind[N_ind[i]]];
            }
            else {
                alpha_i = mu_alpha[condition[i]];
                v_lower_i = mu_v_lower[condition[i]];
                v_upper_i = mu_v_upper[condition[i]];
                tau_i = mu_tau;
            }
            t = (RT[counter]/1000) - tau_i;
            if (t > 0) {
                if (choice[i] == 1) {
                    pdf = race_pdf(t, alpha_i, v_lower_i);
                    cdf = 1 - race_cdf(t, alpha_i, v_upper_i);
                }
                else {
                    pdf = race_pdf(t, alpha_i, v_upper_i);
                    cdf = 1 - race_cdf(t, alpha_i, v_lower_i);
                }
                prob[counter] = pdf*cdf;
                if (prob[counter] < 1e-10) {
                    prob[counter] = 1e-10;
                }
            }
            else {
                prob[counter] = 1e-10;
            }
            counter += 1;

        }
        out = sum(log(prob));
        return out;
    }

    real race_lpdf(real[] RT, int[] choice, int[] N_ind, int[] N_sub_ind, int[] condition, real[] mu_alpha, real[] mu_v_lower, real[] mu_v_upper, real mu_tau, real[,] alpha, real[,] v_lower, real[,] v_upper, real[] tau){

        real t;
        vector[size(RT)] prob;
        real cdf;
        real pdf;
        real out;
        int counter;
        real alpha_i;
        real v_upper_i;
        real v_lower_i;
        real tau_i;

        for (i in 1:size(RT)) {
            if (N_sub_ind[N_ind[i]] != 0) {
                alpha_i = alpha[N_sub_ind[N_ind[i]], condition[i]];
                v_lower_i = v_lower[N_sub_ind[N_ind[i]], condition[i]];
                v_upper_i = v_upper[N_sub_ind[N_ind[i]], condition[i]];
                tau_i = tau[N_sub_ind[N_ind[i]]];
            }
            else {
                alpha_i = mu_alpha[condition[i]];
                v_lower_i = mu_v_lower[condition[i]];
                v_upper_i = mu_v_upper[condition[i]];
                tau_i = mu_tau;
            }
            t = (RT[i]/1000) - tau_i;
            if (t > 0) {
                if (choice[i] == 1) {
                    pdf = race_pdf(t, alpha_i, v_lower_i);
                    cdf = 1 - race_cdf(t, alpha_i, v_upper_i);
                }
                else {
                    pdf = race_pdf(t, alpha_i, v_upper_i);
                    cdf = 1 - race_cdf(t, alpha_i, v_lower_i);
                }
                prob[i] = pdf*cdf;
                if (prob[i] < 1e-10) {
                    prob[i] = 1e-10;
                }
            }
            else {
                prob[i] = 1e-10;
            }
        }
        out = sum(log(prob));
        return out;
    }
}
data {
	int<lower=1> N; // total unique subjects
    int<lower=1> T; // total number of all observations
    int<lower=1> N_sub; // total subjects who get subject-level estimates 
    int N_ind[T]; // index for each unique subject
    int N_sub_ind[N]; // index for each subject who gets subject-level estimates
    int N_cond; // total number of conditions/blocks
    int grainsize; // grainsize
    int condition[T];
    real RT[T];
    int choice[T];
    int sub_id[T];
}
parameters {
    // Hyperparameter means
    real<lower=0> mu_alpha[N_cond];
    real<lower=0.1> mu_tau;
    real<lower=0> mu_v_lower[N_cond];
    real<lower=0> mu_v_upper[N_cond];

    // Hyperparameter variances
    real<lower=0> sigma_alpha[N_cond];
    real<lower=0> sigma_tau;
    real<lower=0> sigma_v_lower[N_cond];
    real<lower=0> sigma_v_upper[N_cond];

    // Individual parameters
    real<lower=0> alpha[N_sub, N_cond];
	real<lower=0.1> tau[N_sub];
	real<lower=0> v_lower[N_sub, N_cond];
	real<lower=0> v_upper[N_sub, N_cond];
}
model {
    // Hyperparameter means
    mu_tau ~ normal(.5, .5)T[0,];
    sigma_tau ~ gamma(1, 1);
    for (i in 1:N_sub) {
        tau[i] ~ normal(mu_tau, sigma_tau)T[0,];
    }
    for (j in 1:N_cond) {
        mu_alpha[j] ~ normal(.5, 1)T[0,];
        mu_v_lower[j] ~ normal(2, 1)T[0,];
        mu_v_upper[j] ~ normal(2, 1)T[0,];
        sigma_alpha[j] ~ gamma(1, 1);
        sigma_v_lower[j] ~ gamma(1, 1);
        sigma_v_upper[j] ~ gamma(1, 1);
        for (k in 1:N_sub) {
            alpha[k, j] ~ normal(mu_alpha[j], sigma_alpha[j])T[0,];
            v_lower[k, j] ~ normal(mu_v_lower[j], sigma_v_lower[j])T[0,];
            v_upper[k, j] ~ normal(mu_v_upper[j], sigma_v_upper[j])T[0,];
        }
    }
    target += reduce_sum(race_reduce, RT, 1, choice, N_ind, N_sub_ind, condition, mu_alpha, mu_v_lower, 
                        mu_v_upper, mu_tau, alpha, v_lower, v_upper, tau);
}
generated quantities {
    real log_lik;
    log_lik = race_lpdf(RT | choice, N_ind, N_sub_ind, condition, mu_alpha, mu_v_lower, 
                        mu_v_upper, mu_tau, alpha, v_lower, v_upper, tau);
}