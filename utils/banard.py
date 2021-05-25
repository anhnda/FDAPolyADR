from dataclasses import make_dataclass
from collections import namedtuple
from dataclasses import make_dataclass
import numpy as np
import warnings
from itertools import combinations
import scipy.stats
from scipy.optimize import shgo
from scipy.special import gamma, kv, gammaln



def _compute_log_combinations(n):
    """Compute all log combination of C(n, k)."""
    gammaln_arr = gammaln(np.arange(n + 1) + 1)
    return gammaln(n + 1) - gammaln_arr - gammaln_arr[::-1]

def _get_binomial_log_p_value_with_nuisance_param(
    nuisance_param, x1_sum_x2, x1_sum_x2_log_comb, index_arr
):
    r"""
    Compute the log pvalue in respect of a nuisance parameter considering
    a 2x2 sample space.
    Parameters
    ----------
    nuisance_param : float
        nuisance parameter used in the computation of the maximisation of
        the p-value. Must be between 0 and 1
    x1_sum_x2 : ndarray
        Sum of x1 and x2 inside barnard_exact
    x1_sum_x2_log_comb : ndarray
        sum of the log combination of x1 and x2
    index_arr : ndarray of boolean
    Returns
    -------
    p_value : float
        Return the maximum p-value considering every nuisance paramater
        between 0 and 1
    Notes
    -----
    Barnard exact test iterate over a nuisance parameter
    :math:`\pi \in [0, 1]` to find the maximum p-value. To search this
    maxima, this function return the negative log pvalue with respect to the
    nuisance parameter passed in params. This negative log p-value is then
    used in `shgo` to find the minimum negative pvalue which is our maximum
    pvalue.
    Also, to compute the different combination used in the
    p-values' computation formula, this function uses `gammaln` which is
    more tolerant for large value than `scipy.special.comb`. `gammaln` gives
    a log combination. For the little precision loss, performances are
    improved a lot.
    """
    t1, t2 = x1_sum_x2.shape
    n = t1 + t2 - 2
    with np.errstate(divide="ignore", invalid="ignore"):
        log_nuisance = np.log(
            nuisance_param,
            out=np.zeros_like(nuisance_param),
            where=nuisance_param >= 0,
        )
        log_1_minus_nuisance = np.log(
            1 - nuisance_param,
            out=np.zeros_like(nuisance_param),
            where=1 - nuisance_param >= 0,
        )

        nuisance_power_x1_x2 = log_nuisance * x1_sum_x2
        nuisance_power_x1_x2[(x1_sum_x2 == 0)[:, :]] = 0

        nuisance_power_n_minus_x1_x2 = log_1_minus_nuisance * (n - x1_sum_x2)
        nuisance_power_n_minus_x1_x2[(x1_sum_x2 == n)[:, :]] = 0

        tmp_log_values_arr = (
            x1_sum_x2_log_comb
            + nuisance_power_x1_x2
            + nuisance_power_n_minus_x1_x2
        )

    tmp_values_from_index = tmp_log_values_arr[index_arr]

    # To avoid dividing by zero in log function and getting inf value,
    # values are centered according to the max
    max_value = tmp_values_from_index.max()

    # To have better result's precision, the log pvalue is taken here.
    # Indeed, pvalue is included inside [0, 1] interval. Passing the
    # pvalue to log makes the interval a lot bigger ([-inf, 0]), and thus
    # help us to achieve better precision
    with np.errstate(divide="ignore", invalid="ignore"):
        log_probs = np.exp(tmp_values_from_index - max_value).sum()
        log_pvalue = max_value + np.log(
            log_probs,
            out=np.full_like(log_probs, -np.inf),
            where=log_probs > 0,
        )

    # Since shgo find the minima, minus log pvalue is returned
    return -log_pvalue



def barnard_exact(table, alternative="two-sided", pooled=True, n=32):
    r"""Perform a Barnard exact test on a 2x2 contingency table.
    Parameters
    ----------
    table : array_like of ints
        A 2x2 contingency table.  Elements should be non-negative integers.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the null and alternative hypotheses. Default is 'two-sided'.
        Please see explanations in the Notes section below.
    pooled : bool, optional
        Whether to compute score statistic with pooled variance (as in
        Student's t-test, for example) or unpooled variance (as in Welch's
        t-test). Default is ``True``.
    n : int, optional
        Number of sampling points used in the construction of the sampling
        method. Note that this argument will automatically be converted to
        the next higher power of 2 since `scipy.stats.qmc.Sobol` is used to
        select sample points. Default is 32. Must be positive. In most cases,
        32 points is enough to reach good precision. More points comes at
        performance cost.
    Returns
    -------
    ber : BarnardExactResult
        A result object with the following attributes.
        statistic : float
            The Wald statistic with pooled or unpooled variance, depending
            on the user choice of `pooled`.
        pvalue : float
            The probability of obtaining a test statistic at least as extreme
            as the one observed under the null hypothesis.
    See Also
    --------
    chi2_contingency : Chi-square test of independence of variables in a
        contingency table.
    fisher_exact : Fisher exact test on a 2x2 contingency table.
    Notes
    -----
    Barnard's test is an exact test used in the analysis of contingency
    tables. It examines the association of two categorical variables, and
    is a more powerful alternative than Fisher's exact test
    for 2x2 contingency tables.
    Let's define :math:`X_0` a 2x2 matrix representing the observed sample,
    where each column stores the binomial experiment, as in the example
    below. Let's also define :math:`p_1, p_2` the theoretical binomial
    probabilities for  :math:`x_{11}` and :math:`x_{12}`. When using
    Barnard exact test, we can assert three different null hypotheses :
    - :math:`H_0 : p_1 \geq p_2` versus :math:`H_1 : p_1 < p_2`,
      with `alternative` = "less"
    - :math:`H_0 : p_1 \leq p_2` versus :math:`H_1 : p_1 > p_2`,
      with `alternative` = "greater"
    - :math:`H_0 : p_1 = p_2` versus :math:`H_1 : p_1 \neq p_2`,
      with `alternative` = "two-sided" (default one)
    In order to compute Barnard's exact test, we are using the Wald
    statistic [3]_ with pooled or unpooled variance.
    Under the default assumption that both variances are equal
    (``pooled = True``), the statistic is computed as:
    .. math::
        T(X) = \frac{
            \hat{p}_1 - \hat{p}_2
        }{
            \sqrt{
                \hat{p}(1 - \hat{p})
                (\frac{1}{c_1} +
                \frac{1}{c_2})
            }
        }
    with :math:`\hat{p}_1, \hat{p}_2` and :math:`\hat{p}` the estimator of
    :math:`p_1, p_2` and :math:`p`, the latter being the combined probability,
    given the assumption that :math:`p_1 = p_2`.
    If this assumption is invalid (``pooled = False``), the statistic is:
    .. math::
        T(X) = \frac{
            \hat{p}_1 - \hat{p}_2
        }{
            \sqrt{
                \frac{\hat{p}_1 (1 - \hat{p}_1)}{c_1} +
                \frac{\hat{p}_2 (1 - \hat{p}_2)}{c_2}
            }
        }
    The p-value is then computed as:
    .. math::
        \sum
            \binom{c_1}{x_{11}}
            \binom{c_2}{x_{12}}
            \pi^{x_{11} + x_{12}}
            (1 - \pi)^{t - x_{11} - x_{12}}
    where the sum is over all  2x2 contingency tables :math:`X` such that:
    * :math:`T(X) \leq T(X_0)` when `alternative` = "less",
    * :math:`T(X) \geq T(X_0)` when `alternative` = "greater", or
    * :math:`T(X) \geq |T(X_0)|` when `alternative` = "two-sided".
    Above, :math:`c_1, c_2` are the sum of the columns 1 and 2,
    and :math:`t` the total (sum of the 4 sample's element).
    The returned p-value is the maximum p-value taken over the nuisance
    parameter :math:`\pi`, where :math:`0 \leq \pi \leq 1`.
    This function's complexity is :math:`O(n c_1 c_2)`, where `n` is the
    number of sample points.
    References
    ----------
    .. [1] Barnard, G. A. "Significance Tests for 2x2 Tables". *Biometrika*.
           34.1/2 (1947): 123-138. :doi:`dpgkg3`
    .. [2] Mehta, Cyrus R., and Pralay Senchaudhuri. "Conditional versus
           unconditional exact tests for comparing two binomials."
           *Cytel Software Corporation* 675 (2003): 1-5.
    .. [3] "Wald Test". *Wikipedia*. https://en.wikipedia.org/wiki/Wald_test
    Examples
    --------
    An example use of Barnard's test is presented in [2]_.
        Consider the following example of a vaccine efficacy study
        (Chan, 1998). In a randomized clinical trial of 30 subjects, 15 were
        inoculated with a recombinant DNA influenza vaccine and the 15 were
        inoculated with a placebo. Twelve of the 15 subjects in the placebo
        group (80%) eventually became infected with influenza whereas for the
        vaccine group, only 7 of the 15 subjects (47%) became infected. The
        data are tabulated as a 2 x 2 table::
                Vaccine  Placebo
            Yes     7        12
            No      8        3
    When working with statistical hypothesis testing, we usually use a
    threshold probability or significance level upon which we decide
    to reject the null hypothesis :math:`H_0`. Suppose we choose the common
    significance level of 5%.
    Our alternative hypothesis is that the vaccine will lower the chance of
    becoming infected with the virus; that is, the probability :math:`p_1` of
    catching the virus with the vaccine will be *less than* the probability
    :math:`p_2` of catching the virus without the vaccine.  Therefore, we call
    `barnard_exact` with the ``alternative="less"`` option:
    >>> import scipy.stats as stats
    >>> res = stats.barnard_exact([[7, 12], [8, 3]], alternative="less")
    >>> res.statistic
    -1.894...
    >>> res.pvalue
    0.03407...
    Under the null hypothesis that the vaccine will not lower the chance of
    becoming infected, the probability of obtaining test results at least as
    extreme as the observed data is approximately 3.4%. Since this p-value is
    less than our chosen significance level, we have evidence to reject
    :math:`H_0` in favor of the alternative.
    Suppose we had used Fisher's exact test instead:
    >>> _, pvalue = stats.fisher_exact([[7, 12], [8, 3]], alternative="less")
    >>> pvalue
    0.0640...
    With the same threshold significance of 5%, we would not have been able
    to reject the null hypothesis in favor of the alternative. As stated in
    [2]_, Barnard's test is uniformly more powerful than Fisher's exact test
    because Barnard's test does not condition on any margin. Fisher's test
    should only be used when both sets of marginals are fixed.
    """
    if n <= 0:
        raise ValueError(
            "Number of points `n` must be strictly positive, "
            f"found {n!r}"
        )

    table = np.asarray(table, dtype=np.int64)

    if not table.shape == (2, 2):
        raise ValueError("The input `table` must be of shape (2, 2).")

    if np.any(table < 0):
        raise ValueError("All values in `table` must be nonnegative.")

    if 0 in table.sum(axis=0):
        # If both values in column are zero, the p-value is 1 and
        # the score's statistic is NaN.
        return np.nan, 1.0

    total_col_1, total_col_2 = table.sum(axis=0)

    x1 = np.arange(total_col_1 + 1, dtype=np.int64).reshape(-1, 1)
    x2 = np.arange(total_col_2 + 1, dtype=np.int64).reshape(1, -1)

    # We need to calculate the wald statistics for each combination of x1 and
    # x2.
    p1, p2 = x1 / total_col_1, x2 / total_col_2

    if pooled:
        p = (x1 + x2) / (total_col_1 + total_col_2)
        variances = p * (1 - p) * (1 / total_col_1 + 1 / total_col_2)
    else:
        variances = p1 * (1 - p1) / total_col_1 + p2 * (1 - p2) / total_col_2

    # To avoid warning when dividing by 0
    with np.errstate(divide="ignore", invalid="ignore"):
        wald_statistic = np.divide((p1 - p2), np.sqrt(variances))

    wald_statistic[p1 == p2] = 0  # Removing NaN values

    wald_stat_obs = wald_statistic[table[0, 0], table[0, 1]]

    if alternative == "two-sided":
        index_arr = np.abs(wald_statistic) >= abs(wald_stat_obs)
    elif alternative == "less":
        index_arr = wald_statistic <= wald_stat_obs
    elif alternative == "greater":
        index_arr = wald_statistic >= wald_stat_obs
    else:
        msg = (
            "`alternative` should be one of {'two-sided', 'less', 'greater'},"
            f" found {alternative!r}"
        )
        raise ValueError(msg)

    x1_sum_x2 = x1 + x2

    x1_log_comb = _compute_log_combinations(total_col_1)
    x2_log_comb = _compute_log_combinations(total_col_2)
    x1_sum_x2_log_comb = x1_log_comb[x1] + x2_log_comb[x2]

    result = shgo(
        _get_binomial_log_p_value_with_nuisance_param,
        args=(x1_sum_x2, x1_sum_x2_log_comb, index_arr),
        bounds=((0, 1),),
        n=n,
        sampling_method="sobol",
    )

    # result.fun is the negative log pvalue and therefore needs to be
    # changed before return
    p_value = np.clip(np.exp(-result.fun), a_min=0, a_max=1)
    return wald_stat_obs, p_value

