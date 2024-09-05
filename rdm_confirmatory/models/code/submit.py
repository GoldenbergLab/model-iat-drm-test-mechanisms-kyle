import pickle
import cmdstanpy
import arviz as az
import warnings
warnings.filterwarnings("ignore")

filename = 'defendattack'
print(filename)

with open('pickles/%s.pkl' % filename, 'rb') as f:
    data = pickle.load(f)

nsubs = data['N']

print(nsubs)

rdm_sm = cmdstanpy.CmdStanModel(
    stan_file='raceddm_new.stan',
    cpp_options={"STAN_THREADS": True},
)

rdm_fit = rdm_sm.sample(data, iter_warmup=1000, iter_sampling=1000, inits=0, chains=4, parallel_chains=4, threads_per_chain=4, seed=101, refresh=1, 
                            save_warmup=1, output_dir='./output/%s' % filename, max_treedepth=10, adapt_delta=0.95)

azdf = az.from_cmdstanpy(posterior=rdm_fit, log_likelihood={'log_lik': 'log_lik'})
az.to_netcdf(azdf, 'raceddm_%s.nc' % filename)
