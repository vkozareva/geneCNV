import logging
import numpy as np
import scipy.optimize

# note that the log-likelihood here is negative since we want to maximize it
# initial reshaping and flattening at end of each function are for compatibility with scipy.optimize.fmin_cg
def _loglike_vcond(z, n, edf, mu, inv_cov):
    """Calculates the negative log likelihood of z conditional on the data and current values of
    mu and cov"""
    z = z.reshape((-1,1))

    neg_loglike = (-(n * np.dot(z.flatten(), edf.flatten()) - n * np.log(1 + np.sum(np.exp(z))) -
                   0.5 *np.dot(np.dot((z - mu).T, inv_cov), (z - mu))))

    return neg_loglike[0][0]

# gradient of the conditional log-likelihood
def _grad_vcond(z, n, edf, mu, inv_cov):
    z = z.reshape((-1,1))
    max_z = np.amax(z)
    pz = np.exp(z - max_z)
    pz = pz / (np.exp(-max_z) + np.sum(pz))

    grad = -(n * (edf - pz) - np.dot(inv_cov,(z - mu)))
    return grad.flatten()

# based on matlab code from https://www.mathworks.com/matlabcentral/fileexchange/11275-hlnfit
def _conditional_mode(y, mu, cov, EMiter, initialz):
    """Calculates the conditional mode of z given data for a single individual and current estimates
    of mu and cov
    y -- kx1 array of read counts in each exon for single individual (histogram of multinomial draws)
    mu -- currrent estimate of mu, k-1 x 1 array
    cov -- current estimate of cov, k-1 x k-1 array
    initialz -- initial estimate for z, k-1 x 1 array

    Some important intermediate computations:
    pz -- the logistically normalized probabilities -- computed by exponentiating z and normalizing by sum

    Returns
    mu_hat -- conditional mode of z, k-1 x 1 array
    cov_hat -- inverse of the negative hessian calculated using conditional mode mu_hat, k-1 x k-1 array
    loglike -- log-likelihood of z and data after convergence, float
    """
    max_iterations = 25
    tol = 1e-9
    k = len(mu)+1
    initialz = initialz.reshape(-1, 1)
    z = initialz
    # get total number of reads/objects drawn from multinomial
    n = np.sum(y)
    # get empirical distribution function (y vector normalized by sum)
    edf = (y / float(n))[:-1].reshape(-1, 1) # must be column vector
    inv_cov = np.linalg.inv(cov)
    grad = np.zeros((k-1, 1))
    neg_hess = np.zeros((k-1, k-1))

    change = tol+1
    iters = 0
    # authors use previously implemented code for Polak-Ribiere conjugate gradient optimization-- instead of recreating, we use
    # the scipy implemented version here
    args = (n, edf, mu, inv_cov)
    z = scipy.optimize.fmin_cg(_loglike_vcond, z, fprime=_grad_vcond, args=args, maxiter=20, disp=0).reshape((-1, 1))
    z_cg = z

    # Newton-raphson after coarse optimization with Polak-Ribiere
    tries = 0
    while iters < max_iterations and change > tol:
        oldz = z

        # note the addition of the final dimension (z_k = 0) -- also subtracting max to prevent underflow
        max_z = np.amax(z)
        pz = np.exp(z - max_z)
        pz = pz / (np.exp(-max_z) + np.sum(pz))

        # gradient with respect to unnormalized intensities (z vector)
        grad = n * (edf - pz) - np.dot(inv_cov, (z - mu))
        # negative hessian (second derivative)
        neg_hess = inv_cov + n * (np.diag(pz.flatten())) - n * (np.dot(pz, pz.T)) # note outer product here
        z = z + np.dot(np.linalg.inv(neg_hess), grad)

        change = np.amax(np.absolute(z - oldz))
        iters += 1

        # try reinitializing with different values if stuck
        if iters > max_iterations - 1:
            z = initialz + 2 * (np.random.rand(k-1, 1) - 0.5 * np.ones((k-1, 1)))
            iters = 0
            tries += 1

        # something is probably wrong if this is actually happening
        if tries > max(EMiter * 5, 20):
            logging.warning('Giving up on Newton Raphson, change is {}'.format(change))
            z = z_cg  # why reset to original starting value in authors' code? shouldn't it be z_cg
            iters = max_iterations

    mu_hat = z
    cov_hat = np.linalg.inv(neg_hess)
    # conditional log likelihood
    loglike = (n * np.dot(z.flatten(), edf.flatten()) - n * np.log(1 + np.sum(np.exp(z))) -
               0.5 *np.dot(np.dot((z - mu).T, inv_cov), (z - mu)))

    return mu_hat, cov_hat, loglike

def hln_EM(Y, max_iterations=25, tol=1e-6, fit_diag_only=False):
    """Returns estimates for mu and cov after EM on subject data.
    Y -- N x k matrix where k is the number of exons and N is the number of subjects
    """
    m = Y.shape[0]
    k = Y.shape[1]
    mu = np.zeros((k-1, 1))
    cov = np.eye(k-1)
    cov_hat = np.zeros((k-1, k-1, m))
    per_sub_loglikes = np.zeros((m, 1))

    change = tol+1
    iters = 0
    mu_hat = np.tile(mu, [1, m])

    while iters < max_iterations and change > tol:
        oldmu = mu
        oldcov = cov

        # for each subject
        for i in range(m):
            mu_hat[:, i:i+1], cov_hat[:, :, i], per_sub_loglikes[i]=_conditional_mode(Y[i, :], mu, cov,
                                                                                      iters, mu_hat[:, i])

        # update based on maximization with multivariate normal approximation
        mu = (1. / m) * np.sum(mu_hat, axis=1).reshape((-1, 1))
        c_mu_hat = mu_hat - np.tile(mu, [1, m])
        cov = (1. / m) * (np.dot(c_mu_hat, c_mu_hat.T) + np.sum(cov_hat, axis=2))
        if fit_diag_only:
            cov = np.diag(np.diagonal(cov))

        change = np.amax([np.amax(np.absolute(mu - oldmu)), np.amax(np.absolute(cov - oldcov))])
        iters += 1
        logging.info('iteration: {}, change: {}'.format(iters, change))

    return mu, cov


def gen_hln_samples(numdocs, numdraws, mu, cov):
    """Generate probability vectors for multinomial distributions from hierarchical logistic normal model
    Not part of the training flow, but useful for testing and experimentation."""
    X = np.exp(np.random.multivariate_normal(mu, cov, numdocs))
    X = np.concatenate((X, np.ones((numdocs, 1))), axis=1)
    X = X / X.sum(axis=1)[:, None]  # add the extra dimension so matrix division returns properly

    Y = np.zeros((numdocs, len(mu)+1))
    for i, x in enumerate(X):
        Y[i] = np.random.multinomial(numdraws, x)

    return Y, X
