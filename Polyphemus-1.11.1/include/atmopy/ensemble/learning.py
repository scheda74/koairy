# Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
#     Author(s): Boris Mauricette, Vivien Mallet, Gilles Stoltz
#
# This file is part of AtmoPy library, a tool for data processing and
# visualization in atmospheric sciences.
#
# AtmoPy is developed in the INRIA - ENPC joint project-team CLIME and in
# the ENPC - EDF R&D joint laboratory CEREA.
#
# AtmoPy is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# AtmoPy is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# For more information, visit the AtmoPy home page:
#     http://cerea.enpc.fr/polyphemus/atmopy.html


import sys, os
sys.path.insert(0,
                os.path.split(os.path.dirname(os.path.abspath(__file__)))[0])
from ensemble_method import *
sys.path.pop(0)


##################################
# EXPONENTIALLY WEIGHTED AVERAGE #
##################################


class ExponentiallyWeightedAverage(EnsembleMethod):
    """
    This class implements the exponentially weighted average algorithm
    (Cesa-Bianchi & Lugosi, 2006, p. 45).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, Nskip = 1, Nlearning = 1, learning_rate =
                 3.e-6, option = "step", verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type learning_rate: float
        @param learning_rate: Learning rate.
        """
        self.learning_rate = learning_rate
        EnsembleMethod.__init__(self, ens,
                                configuration_file = configuration_file,
                                process = process, statistics = statistics,
                                Nskip = Nskip, Nlearning = Nlearning,
                                option = option, verbose = verbose)


    def Init(self):
        self.initial_weight = ones(self.Nsim, 'd') / float(self.Nsim)


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()

        loss = sum((s - o) ** 2, 1)
        weight = previous_weight * exp(-self.learning_rate * loss)
        weight /= weight.sum()

        self.AcquireWeight(weight)


##########################
# EXPONENTIATED GRADIENT #
##########################


class ExponentiatedGradient(EnsembleMethod):
    """
    This class implements the exponentiated gradient algorithm (Cesa-Bianchi,
    1999; Cesa-Bianchi & Lugosi, 2006, p. 23).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, extended = False, U = 1., Nskip = 1,
                 Nlearning = 1, learning_rate = 2.e-5, option = "step",
                 verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type learning_rate: float
        @param learning_rate: Learning rate.
        """
        self.learning_rate = learning_rate
        EnsembleMethod.__init__(self, ens,
                                configuration_file = configuration_file,
                                process = process, statistics = statistics,
                                Nskip = Nskip, Nlearning = Nlearning,
                                option = option, extended = extended,
                                verbose = verbose, U = U)


    def Init(self):
        self.initial_weight = ones((self.Nsim), 'd') / float(self.Nsim)


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()

        loss = 2. * sum((dot(previous_weight, s) - o) * s, 1)
        weight = previous_weight * exp(-self.learning_rate * loss)
        weight /= weight.sum()

        self.AcquireWeight(weight)


######################################
# EXPONENTIATED GRADIENT WITH WINDOW #
######################################


class ExponentiatedGradientWindow(ExponentiatedGradient):
    """
    This class implements a modified exponentiated gradient algorithm (Mallet,
    Mauricette, and Stoltz, 2007).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, extended = False, U = 1., Nskip = 1,
                 Nlearning = 1, learning_rate = 2.e-5, Nkeep = 20,
                 option = "step", verbose = False):
        """
        See documentation of 'ExponentiatedGradient.__init__' for explanations
        about arguments.

        @type Nkeep: integer
        @param Nkeep: number of training steps.
        """
        self.Nkeep = Nkeep
        ExponentiatedGradient.__init__(self, ens, configuration_file =
                                       configuration_file, process = process,
                                       statistics = statistics, Nskip = Nskip,
                                       Nlearning = Nlearning, extended =
                                       extended, U = U, option = option,
                                       learning_rate = learning_rate, verbose
                                       = verbose)


    def Init(self):
        self.initial_weight = ones((self.Nsim), 'd') / float(self.Nsim)
        self.kept_weight = self.InitialList([])


    def GetTools(self):
        if self.ens.config.concentrations == "peak":
            return self.kept_weight[0]
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            return self.kept_weight[hour]


    def UpdateTools(self, kept_weight):
        if self.ens.config.concentrations == "peak":
            self.kept_weight[0] = kept_weight
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            self.kept_weight[hour] = kept_weight


    def Upkeep(self, kept_weight, previous_weight, loss):
        if len(kept_weight) < self.Nkeep:
            kept_weight.append(exp(-self.learning_rate * loss))
            weight = previous_weight * exp(-self.learning_rate * loss)
        elif len(kept_weight) == self.Nkeep:
            kept_weight.append(exp(-self.learning_rate * loss))
            old_confidence = kept_weight.pop(0)
            weight = previous_weight * exp(-self.learning_rate * loss) \
                     / old_confidence
        weight /= weight.sum()
        return weight, kept_weight


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()
        kept_weight = self.GetTools()

        loss = 2. * sum((dot(previous_weight, s) - o) * s, 1)
        weight, kept_weight = self.Upkeep(kept_weight, previous_weight, loss)

        self.UpdateTools(kept_weight)
        self.AcquireWeight(weight)


#####################################
# EXPONENTIATED GRADIENT DISCOUNTED #
#####################################


class ExponentiatedGradientDiscounted(EnsembleMethod):
    """
    This class implements the exponentiated gradient algorithm with discount
    losses (Mallet, Mauricette, and Stoltz, 2007).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, extended = False, U = 1., Nskip = 1,
                 Nlearning = 1, learning_rate = 1.2e-4, forget_rate = 1.,
                 p1 = 0.5, p2 = 1., option = "step", verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type learning_rate: float
        @param learning_rate: Learning rate.
        @type forget_rate: float
        @param forget_rate: Forget rate of past confidence.
        @type p1: float
        @param p1: Power of the learning rate.
        @type p2: float
        @param p2: Power of the forget rate.
        """
        self.learning_rate = learning_rate
        self.forget_rate = forget_rate
        self.p1 = p1
        self.p2 = p2
        EnsembleMethod.__init__(self, ens,
                                configuration_file = configuration_file,
                                process = process, statistics = statistics,
                                Nskip = Nskip, Nlearning = Nlearning,
                                option = option, extended = extended,
                                verbose = verbose, U = U)


    def Init(self):
        self.initial_weight = ones(self.Nsim, 'd') / float(self.Nsim)
        self.confidence = self.InitialList([])


    def GetTools(self):
        if self.ens.config.concentrations == "peak":
            return self.confidence[0]
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            return self.confidence[hour]


    def UpdateTools(self, confidence):
        if self.ens.config.concentrations == "peak":
            self.confidence[0] = confidence
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            self.confidence[hour] = confidence


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()
        confidence = self.GetTools()

        loss = 2. * sum((dot(previous_weight, s) - o) * s, 1)
        confidence.append(loss)
        T = len(confidence)
        weight = ones(len(previous_weight), 'd')
        for i in range(T):
            weight *= exp(-self.learning_rate * confidence[i]
                           / T ** self.p1 * (1 + self.forget_rate
                                             / (T - i) ** self.p2))
            weight /= weight.sum()

        self.UpdateTools(confidence)
        self.AcquireWeight(weight)


###################################
# ADAPTIVE EXPONENTIATED GRADIENT #
###################################


class ExponentiatedGradientAdaptive(ExponentiatedGradientDiscounted):
    """
    This class implements the exponentiated gradient algorithm
    with adaptive learning rate (Stoltz, Cesa-Bianchi, and Mansour, 2007).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, extended = False, U = 1., Nskip = 1,
                 Nlearning = 1, a = 100., b = 1.,
                 option = "step", verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type a: float
        @param a: A parameter to compute the learning rate.
        @type b: float
        @param b: Another parameter to compute the learning rate.
        """
        self.a = a
        self.b = b
        ExponentiatedGradientDiscounted.__init__(self, ens, configuration_file
                                                 = configuration_file, process
                                                 = process, statistics =
                                                 statistics, Nskip = Nskip,
                                                 Nlearning = Nlearning,
                                                 extended = extended, option =
                                                 option, verbose = verbose, U
                                                 = U)


    def Init(self):
        self.factor = sqrt(2. * sqrt(2.) - 1.) / (exp(1.) - 2.)
        self.initial_weight = ones(self.Nsim,'d') / float(self.Nsim)
        self.confidence = self.InitialList([])
        self.variance = self.InitialList(0.)
        self.bound = self.InitialList(0.)
        self.learning_rate = self.InitialList(1.)


    def GetTools(self):
        if self.ens.config.concentrations == "peak":
            return self.confidence[0], self.variance[0], self.bound[0], \
                   self.learning_rate[0]
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            return self.confidence[hour], self.variance[hour], \
                   self.bound[hour], self.learning_rate[hour]


    def UpdateTools(self, confidence, variance, bound, learning_rate):
        if self.ens.config.concentrations == "peak":
            self.confidence[0] = confidence
            self.variance[0] = variance
            self.bound[0] = bound
            self.learning_rate[0] = learning_rate
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            self.confidence[hour] = confidence
            self.variance[hour] = variance
            self.bound[hour] = bound
            self.learning_rate[hour] = learning_rate


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()
        confidence, variance, bound, learning_rate = self.GetTools()
        loss = 2. * sum((dot(previous_weight, s) - o) * s, 1)
        mean_loss = inner(previous_weight, loss)

        # Adds the variance of the gradient-loss term.
        variance += inner(previous_weight, (loss - mean_loss) ** 2)
        bound = max(bound, abs(loss).max())
        learning_rate = min(self.a / bound,
                            self.factor * self.b
                            * sqrt(log(self.Nsim) / variance))
        confidence.append(loss)
        T = len(confidence)
        weight = ones(len(previous_weight), 'd')
        for i in range(T):
            weight *= exp(-learning_rate * confidence[i])
            weight /= weight.sum()

        self.UpdateTools(confidence, variance, bound, learning_rate)
        self.AcquireWeight(weight)


########
# PROD #
########


class Prod(ExponentiatedGradient):
    """
    This class implements the algorithm Prod (Cesa-Bianchi, Mansour, and
    Stoltz, 2007).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, extended = False, U = 1., Nskip = 1,
                 Nlearning = 1, learning_rate = 5.5e-7,
                 option = "step", verbose = False):

        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type learning_rate: float
        @param learning_rate: Learning rate.
        """
        ExponentiatedGradient.__init__(self, ens, configuration_file =
                                       configuration_file, process = process,
                                       statistics = statistics, Nskip = Nskip,
                                       Nlearning = Nlearning, extended =
                                       extended, U = U, option = option,
                                       learning_rate = learning_rate, verbose
                                       = verbose)


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()

        loss = 2. * sum((dot(previous_weight, s) - o) * s, 1)
        weight = previous_weight * (1. - self.learning_rate * loss)
        if (weight < 0.).any():
            raise Exception, "Too large parameter eta"
        weight /= weight.sum()

        self.AcquireWeight(weight)


####################
# GRADIENT DESCENT #
####################


class GradientDescent(EnsembleMethod):
    """
    This class implements algorithm GD (Cesa-Bianchi, 1999).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, Nskip = 1, Nlearning = 1,
                 learning_rate = 4.5e-9, lamb = 1.,
                 option = "step", verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type learning_rate: float
        @param learning_rate: Learning rate.
        @type lamb: float
        @param lamb: Lambda coefficient of initial vector.
        """
        self.lamb = lamb
        self.learning_rate = learning_rate
        EnsembleMethod.__init__(self, ens,
                                configuration_file = configuration_file,
                                process = process, statistics = statistics,
                                Nskip = Nskip, Nlearning = Nlearning,
                                option = option, verbose = verbose)


    def Init(self):
        self.initial_weight = self.lamb * ones((self.Nsim), 'd') \
                              / float(self.Nsim)


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()

        loss = 2. * sum((dot(previous_weight, s) - o) * s, 1)
        weight = previous_weight - self.learning_rate * loss

        self.AcquireWeight(weight)


#####################
# GREEDY PROJECTION #
#####################


def projection_simplex(v):
    """
    Projection of the vector onto the simplex of probability distributions.

    @type v: array
    @param v: array.
    """
    v = -v
    li = argsort(v)
    v.sort()
    v = -v
    vc = 0.
    m = 1
    while vc <= 1. and m < len(v):
        vc += m * (v[m-1] - v[m])
        m += 1
    m -= 1
    z = v[:m]
    la = 1. / m * (1 - sum(z))
    z = z + la
    z = concatenate([z, zeros(len(v) - m)], 0)
    Inv = [list(li).index(j) for j in range(len(v))]
    return z[Inv]


def projection_cubic(v, r):
    """
    Projection of the vector onto the cube of radius r.

    @type v: array
    @param v: array.
    @type r: float
    @param r: radius of the cube.
    """
    return minimum(maximum(-r, v), r)


def projection_L2(v, r):
    """
    Projection of the vector onto the ball of radius r.

    @type v: array
    @param v: array.
    @type r: float
    @param r: radius of the ball.
    """
    Norm = sqrt((v * v).sum())
    if Norm > r:
        return r * v / Norm
    else:
        return v


class Zink(EnsembleMethod):
    """
    This class implements the greedy projection gradient algorithm (Zinkevich,
    2003).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, extended = False, U = 1., Nskip = 1,
                 Nlearning = 1, learning_rate = 1.e-6, radius = 1.,
                 projection_function = projection_simplex,
                 option = "step", verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type learning_rate: float
        @param learning_rate: Learning rate.
        @type projection_function: function
        @param projection_function: Projection on simplex or ball
        @type radius: float
        @param radius: radius of the ball where to project.
        """
        self.learning_rate = learning_rate
        self.radius = radius
        self.projection_function = projection_function
        EnsembleMethod.__init__(self, ens,
                                configuration_file = configuration_file,
                                process = process, statistics = statistics,
                                Nskip = Nskip, Nlearning = Nlearning,
                                option = option, extended = extended,
                                verbose = verbose, U = U)


    def Init(self):
        self.initial_weight = ones((self.Nsim), 'd') / float(self.Nsim)
        if not(self.projection_function == projection_simplex) \
               and self.extended:
            raise Exception, \
                  "'extended' option only allowed with projection_simplex."


    def UpdateWeight(self, s, o):
        import inspect
        previous_weight = self.GetPreviousWeight()
        loss = 2. * sum((dot(previous_weight, s) - o) * s, 1)

        if len(inspect.getargspec(self.projection_function)[0]) == 1:
            weight = self.projection_function(previous_weight
                                              - self.learning_rate * loss)
        else:
            weight = self.projection_function(previous_weight
                                              - self.learning_rate * loss,
                                              self.radius)

        self.AcquireWeight(weight)


####################
# RIDGE REGRESSION #
####################


class RidgeRegression(EnsembleMethod):
    """
    This class implements a modified ridge-regression algorithm
    (Cesa-Bianchi & Lugosi, 2006, p. 317).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, Nskip = 1, Nlearning = 1,
                 penalization = 1., option = "step", verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type penalization: float
        @param penalization: penalization of the norm of the weights.
        """
        self.penalization = penalization
        EnsembleMethod.__init__(self, ens,
                                configuration_file = configuration_file,
                                process = process, statistics = statistics,
                                Nskip = Nskip, Nlearning = Nlearning,
                                option = option, verbose = verbose)


    def Init(self):
        self.initial_weight = zeros(self.Nsim)
        self.A = self.InitialList(self.penalization * identity(self.Nsim))


    def GetTools(self):
        if self.ens.config.concentrations == "peak":
            return self.A[0]
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            return self.A[hour]


    def UpdateTools(self, A):
        if self.ens.config.concentrations == "peak":
            self.A[0] = A
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            self.A[hour] = A


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()
        A = self.GetTools()

        Nobs = len(o)
        bt = zeros(self.Nsim)
        for i in range(Nobs):
            x_it = s[:, i]
            A += outer(x_it, x_it)
            bt += (inner(x_it, previous_weight) - o[i]) * x_it
        # Warning: A grows quickly (but linearly).
        weight = previous_weight - dot(scipy.linalg.pinv2(A), bt)

        self.UpdateTools(A)
        self.AcquireWeight(weight)


#############################
# RIDGE REGRESSION MODIFIED #
#############################


class RidgeRegressionModified(RidgeRegression):
    """
    This class implements a modified ridge-regression algorithm
    (Cesa-Bianchi & Lugosi, 2006, chap 11.8).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, Nskip = 0, Nlearning = 0,
                 penalization = 1.e6, option = "step", verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type penalization: float
        @param penalization: penalization of the norm of the weights.
        """
        RidgeRegression.__init__(self, ens,
                                 configuration_file = configuration_file,
                                 process = process, statistics = statistics,
                                 Nskip = Nskip, Nlearning = Nlearning,
                                 option = option, penalization = penalization,
                                 verbose = verbose)


    def Init(self):
        self.initial_weight = zeros(self.Nsim, 'd')
        self.A = self.InitialList(self.penalization * identity(self.Nsim))
        self.b = self.InitialList(zeros(self.Nsim, 'd'))


    def GetTools(self):
        if self.ens.config.concentrations == "peak":
            return self.A[0], self.b[0]
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            return self.A[hour], self.b[hour]


    def UpdateTools(self, A, b):
        if self.ens.config.concentrations == "peak":
            self.A[0] = A
            self.b[0] = b
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            self.A[hour] = A
            self.b[hour] = b


    def UpdateWeight(self, s, o):
        prev_weight = self.GetPreviousWeight()
        A, b = self.GetTools()

        index = range(len(o))
        tmp_mat = [outer(s[:, i], s[:, i]) for i in index]
        A += reduce(add, tmp_mat)

        weight = dot(scipy.linalg.pinv2(A), b)
        # 'b' is updated after forecast.
        b += reduce(add, map(multiply,
                             [s[:,i] for i in index], o))

        self.UpdateTools(A, b)
        self.AcquireWeight(weight)


###############################
# RIDGE REGRESSION DISCOUNTED #
###############################


class RidgeRegressionDiscounted(RidgeRegressionModified):
    """
    This class implements a modified ridge-regression algorithm
    (Cesa-Bianchi & Lugosi, 2006, p. 317).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, Nskip = 1, Nlearning = 1,
                 option = "step", penalization = 1000.,
                 forget_rate = 100., p1 = 2., verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type penalization: float
        @param penalization: penalization of the norm of the weights.
        @type forget_rate: float
        @param forget_rate: Forget rate of past confidence.
        @type p1: float
        @param p1: power in the forget coefficient.
        """
        self.forget_rate = forget_rate
        self.p1 = p1
        RidgeRegressionModified.__init__(self, ens, configuration_file =
                                         configuration_file, process =
                                         process, statistics = statistics,
                                         Nskip = Nskip, Nlearning = Nlearning,
                                         penalization = penalization, option =
                                         option, verbose = verbose)


    def Init(self):
        self.initial_weight = zeros(self.Nsim)
        self.A = self.InitialList([zeros([self.Nsim, self.Nsim])])
        self.b = self.InitialList([zeros(self.Nsim)])


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()
        At, b = self.GetTools()

        tmp_mat = [outer(s[:, i], s[:, i]) for i in range(len(o))]
        tmp_b = reduce(add, map(multiply, o, [s[:,i] for i in range(len(o))]))
        At.append(reduce(add, tmp_mat))
        b.append(tmp_b)

        T = len(At)
        forget = [(self.forget_rate / ((T - i) ** self.p1) + 1.)
                  for i in range(T)]
        A = reduce(add, map(multiply, At, forget))
        bt = reduce(add, map(multiply, b, forget))
        A += self.penalization * identity(self.Nsim)
        weight = dot(scipy.linalg.pinv2(A), bt)

        self.UpdateTools(At, b)
        self.AcquireWeight(weight)


###########################
# RIDGE REGRESSION WINDOW #
###########################


class RidgeRegressionWindow(RidgeRegressionModified):
    """
    This class implements a modified ridge-regression algorithm
    (Cesa-Bianchi & Lugosi, 2006, p. 317).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, Nskip = 1, Nlearning = 1,
                 penalization = 100., Nkeep = 45, option = "step",
                 verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type penalization: float
        @param penalization: penalization of the norm of the weights.
        """
        self.Nkeep = Nkeep
        RidgeRegressionModified.__init__(self, ens, configuration_file =
                                         configuration_file, process =
                                         process, statistics = statistics,
                                         Nskip = Nskip, Nlearning = Nlearning,
                                         option = option, verbose = verbose,
                                         penalization = penalization)


    def Init(self):
        self.initial_weight = zeros(self.Nsim)
        self.A = self.InitialList([])
        self.b = self.InitialList([])


    def Upkeep(self, A, At, b, bt):
        if len(A) < self.Nkeep:
            A.append(At.copy())
            b.append(bt.copy())
        elif len(A) == self.Nkeep:
            A.append(At.copy())
            b.append(bt.copy())
            del A[0]
            del b[0]
        else:
            raise Exception, "Length of 'A' strictly greater than 'Nkeep'!"
        return A, b


    def UpdateWeight(self, s, o):
        period = self.GetLearningDates()
        s, o = self.CollectData(period)
        A, b = self.GetTools()

        Nobs = len(o)
        At = zeros([self.Nsim, self.Nsim], 'd')
        bt = zeros(self.Nsim, 'd')
        for i in range(Nobs):
            x_it = s[:, i]
            At += outer(x_it, x_it)
            bt += o[i] * x_it

        A, b = self.Upkeep(A, At, b, bt)
        Atmp = reduce(add, A)
        btmp = reduce(add, b)
        Atmp += self.penalization * identity(self.Nsim)
        weight = dot(scipy.linalg.pinv2(Atmp), btmp)

        self.UpdateTools(A, b)
        self.AcquireWeight(weight)


###########
# MIXTURE #
###########


def uniform_simplex(N):
    """
    Simulation of a uniform probability on the simplex.

    @type N: integer
    @param N: The number of dimensions.
    @rtype: 1D array
    @return: An element of the simplex, randomly generated.
    """
    P = []
    list = [0., 1.]
    for i in range(N - 1):
        list.append(random.uniform(0, 1))
    list = sort(list)
    for i in range(N):
        P.append(list[i+1] - list[i])
    return array(P)


class Mixture(EnsembleMethod):
    """
    This class implements the exponentially weighted average mixture algorithm
    (Cesa-Bianchi & Lugosi, 2006, p. 48).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, extended = False, U = 1., Nskip = 1,
                 Nlearning = 1, learning_rate = 1.e-4, Napprox = 5000,
                 option = "step", verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type learning_rate: float
        @param learning_rate: Learning rate.
        @type Napprox: integer
        @param Napprox: Number of points to approximate the integral over the
        simplex with a trivial quadrature formula.
        """
        self.Napprox = Napprox
        self.learning_rate = learning_rate
        EnsembleMethod.__init__(self, ens,
                                configuration_file = configuration_file,
                                process = process, statistics = statistics,
                                Nskip = Nskip, Nlearning = Nlearning,
                                option = option, extended = extended,
                                verbose = verbose, U = U)


    def Init(self):
        self.initial_weight = ones(self.Nsim, 'd') / float(self.Nsim)
        self.confidence = self.InitialList(ones(self.Napprox, 'd')
                                           / float(self.Napprox))
        if self.extended:
            self.interpolweights = zeros([2 * self.Nsim, self.Napprox])
            # Keeps the points of quadrature in the simplex.
            for i in range(self.Napprox):
                self.interpolweights[:, i] = uniform_simplex(2 * self.Nsim)
        else:
            self.interpolweights = zeros([self.Nsim, self.Napprox])
            # Keeps the points of quadrature in the simplex.
            for i in range(self.Napprox):
                self.interpolweights[:, i] = uniform_simplex(self.Nsim)


    def GetTools(self):
        if self.ens.config.concentrations == "peak":
            return self.confidence[0]
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            return self.confidence[hour]


    def UpdateTools(self, confidence):
        if self.ens.config.concentrations == "peak":
            self.confidence[0] = confidence
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            self.confidence[hour] = confidence


    def UpdateWeight(self, s, o):
        confidence = self.GetTools()

        for i in range(self.Napprox):
            loss = exp(-self.learning_rate
                       * sum((dot(self.interpolweights[:, i], s) - o) ** 2))
            confidence[i] *= loss
        confidence /= confidence.sum()
        weight = sum(confidence * self.interpolweights, 1)

        self.UpdateTools(confidence)
        self.AcquireWeight(weight)


###################
# POLYNOMIAL LOSS #
###################


class Polynomial(EnsembleMethod):
    """
    This class implements the polynomially weighted average forecaster
    (Cesa-Bianchi & Lugosi, 2006, p. 12).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, Nskip = 1, Nlearning = 1, power = 2.,
                 option = "step", verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type power: float
        @param power: Power of the potential. It must be strictly greater than
        one.
        """
        self.power = power
        EnsembleMethod.__init__(self, ens, configuration_file =
                                configuration_file, process = process,
                                statistics = statistics, Nskip = Nskip,
                                Nlearning = Nlearning, option = option,
                                verbose = verbose)


    def Init(self):
        self.initial_weight = ones(self.Nsim, 'd') / float(self.Nsim)
        self.regret = self.InitialList(zeros(self.Nsim, 'd'))


    def GetTools(self):
        if self.ens.config.concentrations == "peak":
            return self.regret[0]
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            return self.regret[hour]


    def UpdateTools(self, regret):
        if self.ens.config.concentrations == "peak":
            self.regret[0] = regret
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            self.regret[hour] = regret


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()
        regret = self.GetTools()

        loss = sum((dot(previous_weight, s) - o) ** 2)
        regret += loss - sum((s - o) ** 2, 1)
        # Keeps the positive part of regret.
        weight = ((regret + abs(regret)) / 2.) ** (self.power - 1)
        if (weight == 0).all():
            weight = previous_weight
        else:
            weight /= weight.sum()

        self.UpdateTools(regret)
        self.AcquireWeight(weight)


############################
# GRADIENT POLYNOMIAL LOSS #
############################


class PolynomialGradient(Polynomial):
    """
    This class implements polynomially weighted average forecaster using the
    gradient of the loss (Cesa-Bianchi & Lugosi, 2006, pp. 12 & 23).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, extended = False, U = 1., Nskip = 1,
                 Nlearning = 1, power = 14., option = "step",
                 verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type power: float
        @param power: Power of the potential. It must be strictly greater than
        one.
        """
        self.power = power
        EnsembleMethod.__init__(self, ens, configuration_file =
                                configuration_file, process = process,
                                statistics = statistics, Nskip = Nskip,
                                Nlearning = Nlearning, option = option,
                                extended = extended,
                                verbose = verbose, U = U)


    def Init(self):
        self.initial_weight = ones((self.Nsim), 'd') / float(self.Nsim)
        if self.extended:
            self.regret = self.InitialList(zeros(2 * self.Nsim, 'd'))
        else:
            self.regret = self.InitialList(zeros(self.Nsim, 'd'))


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()
        regret = self.GetTools()

        forecast = dot(previous_weight, s)
        gradloss = 2. * (forecast - o)
        regret += dot(forecast - s, transpose(gradloss))
        weight = ((regret + abs(regret)) / 2) ** (self.power - 1)
        if (weight == 0).all():
            weight = previous_weight
        else:
            weight /= weight.sum()

        self.UpdateTools(regret)
        self.AcquireWeight(weight)


###############
# FIXED-SHARE #
###############


class FixedShare(EnsembleMethod):
    """
    This class implements the fixed-share algorithm
    (Cesa-Bianchi & Lugosi, 2006, Section 5.2).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, Nskip = 1, Nlearning = 1,
                 learning_rate = 1.5e-5, shift = 5.e-2,
                 option = "step", verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type learning_rate: float
        @param learning_rate: Learning rate.
        @type shift: float
        @param shift: Probability to change the leading method.
        """
        self.learning_rate = learning_rate
        self.shift = shift
        EnsembleMethod.__init__(self, ens,
                                configuration_file = configuration_file,
                                process = process, statistics = statistics,
                                Nskip = Nskip, Nlearning = Nlearning,
                                option = option, verbose = verbose)


    def Init(self):
        self.initial_weight = ones(self.Nsim,'d') / float(self.Nsim)
        self.confidence = self.InitialList(ones(self.Nsim, 'd')
                                           / float(self.Nsim))


    def GetTools(self):
        if self.ens.config.concentrations == "peak":
            return self.confidence[0]
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            return self.confidence[hour]


    def UpdateTools(self, confidence):
        if self.ens.config.concentrations == "peak":
            self.confidence[0] = confidence
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            self.confidence[hour] = confidence


    def UpdateWeight(self, s, o):
        confidence = self.GetTools()

        loss = sum((s - o) ** 2, 1)
        weight_tmp = confidence * exp(-self.learning_rate * loss)
        confidence = self.shift * sum(weight_tmp) / float(self.Nsim) \
                     + (1. - self.shift) * weight_tmp
        weight = confidence / confidence.sum()

        self.UpdateTools(confidence)
        self.AcquireWeight(weight)


########################
# FIXED-SHARE GRADIENT #
########################


class FixedShareGradient(FixedShare):
    """
    This class implements the fixed-share algorithm using the gradient of the
    loss (Cesa-Bianchi & Lugosi, 2006, Section 5.2).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, extended = False, U = 1., Nskip = 1,
                 Nlearning = 1, learning_rate = 2.5e-5, shift = 2.e-2,
                 option = "step", verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type learning_rate: float
        @param learning_rate: Learning rate.
        @type shift: float
        @param shift: Probability to change the leading method.
        """
        self.learning_rate = learning_rate
        self.shift = shift
        EnsembleMethod.__init__(self, ens,
                                configuration_file = configuration_file,
                                process = process, statistics = statistics,
                                Nskip = Nskip, Nlearning = Nlearning,
                                option = option, extended = extended,
                                U = U, verbose = verbose)


    def Init(self):
        self.initial_weight = ones(self.Nsim,'d') / float(self.Nsim)
        # Advice: modify the initial confidence to an a priori
        # measure (p. 102).
        if self.extended:
            self.confidence = self.InitialList(ones(2 * self.Nsim, 'd')
                                               / float(2 * self.Nsim))
        else:
            self.confidence = self.InitialList(ones(self.Nsim, 'd')
                                               / float(self.Nsim))


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()
        confidence = self.GetTools()

        loss = 2. * sum((dot(previous_weight, s) - o) * s, 1)
        weight_tmp = confidence * exp(-self.learning_rate * loss)
        confidence = self.shift * sum(weight_tmp) / float(self.Nsim) \
                     + (1. - self.shift) * weight_tmp
        weight = confidence / confidence.sum()

        self.UpdateTools(confidence)
        self.AcquireWeight(weight)


###########################
# VARIABLE-SHARE GRADIENT #
###########################


class VariableShareGradient(FixedShareGradient):
    """
    This class implements the variable-share algorithm using the
    gradient of the loss (Cesa-Bianchi & Lugosi, 2006, Section 5.2).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, extended = False, U = 1., Nskip = 1,
                 Nlearning = 1, learning_rate = 1., shift = 5.e-2,
                 option = "step", verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type learning_rate: float
        @param learning_rate: Learning rate.
        @type shift: float
        @param shift: Probability to change the leading method.
        """
        FixedShareGradient.__init__(self, ens,
                                    configuration_file = configuration_file,
                                    process = process,
                                    statistics = statistics,
                                    Nskip = Nskip,
                                    Nlearning = Nlearning,
                                    extended = extended,
                                    learning_rate = learning_rate,
                                    shift = shift, option = option, U = U,
                                    verbose = verbose)


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()
        confidence = self.GetTools()

        # Loss normalized by L.
        loss = 2. * sum((dot(previous_weight, s) - o) * s, 1) / 1.e6
        # N : number of models, even in extended mode.
        N = s.shape[0]
        weight_tmp = confidence * exp(-self.learning_rate * loss)
        A = zeros([N, N] , 'd')
        for j in range(N):
            for i in range(N):
                if i == j:
                    A[i, j] = (1. - self.shift) ** loss[j]
                else:
                    A[i, j] = (1. - (1. - self.shift) ** loss[j]) \
                              / float(N - 1)
        confidence = dot(A, weight_tmp)
        weight = confidence / confidence.sum()

        self.UpdateTools(confidence)
        self.AcquireWeight(weight)


######################
# ONLINE NEWTON STEP #
######################


def projASimplex(A, y):
    """
    Projection on the simplex under A^-1 distance.

    @type A : array
    @param A: array.
    @type y : array
    @param y: vector to project.
    """
    s = array(zeros(len(y)),  'f')
    x = array(y)
    # Precision.
    prec =  1.e-4
    while (sum(s - x > prec) >  0) or (sum(x - s > prec) >  0):
       x = s
       grad = dot(A, s - y)
       rho = sum(grad * grad) / dot(transpose(grad), dot(A, grad))
       uu = s - rho * grad
       s = projection_simplex(uu)
    return s


class OnlineNewtonStep(EnsembleMethod):
    """
    This class implements the ONS algorithm (see Hazan, Kalai, Kale, and
    Agarwal, 2006).
    """


    def __init__(self, ens, configuration_file = None,
                 process = True, verbose = False,
                 statistics = True, Nskip = 1, Nlearning = 1,
                 beta = 6.e-7, option = "step", mix = 0.2,
                 extended = False, U = 1.):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type beta: float
        @param beta: A parameter in the body of the algorithm.
        @type mix: float
        @param mix: A parameter in [0, 1] to mix the weight with the uniform
        probability.
        """
        self.beta = beta
        self.mix = mix
        EnsembleMethod.__init__(self, ens,
                                configuration_file = configuration_file,
                                process = process, statistics = statistics,
                                Nskip = Nskip, Nlearning = Nlearning,
                                option = option, verbose = verbose, extended =
                                extended, U = U)


    def Init(self):
        self.initial_weight = ones(self.Nsim, 'd') / float(self.Nsim)
        if self.extended:
            N = 2 * self.Nsim
        else:
            N = self.Nsim
        self.A = self.InitialList(zeros([N, N], 'd'))
        self.b = self.InitialList(zeros(N, 'd'))
        self.uniform = ones(N, 'd') / float(N)


    def GetTools(self):
        if self.ens.config.concentrations == "peak":
            return self.A[0], self.b[0]
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            return self.A[hour], self.b[hour]


    def UpdateTools(self, At, bt):
        if self.ens.config.concentrations == "peak":
            self.A[0] = At
            self.b[0] = bt
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            self.A[hour] = At
            self.b[hour] = bt


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()
        At, bt = self.GetTools()

        forecast = dot(previous_weight, s)
        loss = 2. * sum((forecast - o) * s, 1)
        At += outer(loss, loss)
        bt += dot(outer(loss, loss), previous_weight) \
                  - transpose(loss) / self.beta
        weight = (1. - self.mix) \
                 * projASimplex(At, dot(scipy.linalg.pinv2(At), bt)) \
                 + self.mix * self.uniform

        self.UpdateTools(At, bt)
        self.AcquireWeight(weight)


###################
# INTERNAL METHOD #
###################


def i2j(N,i,j):
    """
    Returns the matrix that, applied to a vector, adds the i-th element to the
    j-th element, and sets the i-th element to 0.

    @type N: int
    @param N: The size of the matrix.
    @type i: int
    @param i: The index of the element to be added and to be set to zero.
    @type j: int
    @param j: The index of the element to be augmented.

    @rtype: 2D array
    @return: The proper matrix.
    """
    M = identity(N, dtype = 'd')
    M[i, i] = 0.
    M[j, i] = 1.
    return M


def fixed_point2(Q, prec, w):
    """
    @type Q: square array
    @type prec: float
    @type w: array (vector)
    """
    def MatZero(Q):
        N = Q.shape[0]
        M = transpose(Q) - diag(sum(Q,1))
        return M

    n = Q.shape[0]
    M = MatZero(Q)
    Id = identity(n, dtype = 'd')
    R = M + Id
    weight = w
    S = 1.
    while (S > prec):
        wN = dot(R, weight)
        MC = dot(M, wN)
        # Quick multiplication.
        R = dot(R, R)
        S = max(abs(max(MC)), abs(min(MC)))
        weight = wN
    return weight


class InternalMethod(EnsembleMethod):
    """
    This class implements the template for internal computation of weights
    (Cesa-Bianchi, 1999; Cesa-Bianchi & Lugosi, 2006, p. 23).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, extended = False, U = 1., Nskip = 1,
                 Nlearning = 1, option = "step", precision = 1.e-10, verbose =
                 False, fixed_point_method = fixed_point2):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.
        """
        self.precision = precision
        self.fixed_point_method = fixed_point_method
        EnsembleMethod.__init__(self, ens,
                                configuration_file = configuration_file,
                                process = process, statistics = statistics,
                                Nskip = Nskip, Nlearning = Nlearning,
                                option = option, extended = extended,
                                U = U, verbose = verbose)


    def Init(self):
        self.initial_weight = ones((self.Nsim), 'd') / float(self.Nsim)
        if self.extended:
            self.Q = self.InitialList(ones([2 *self.Nsim , 2 * self.Nsim],
                                           'd')
                                      - identity(2 * self.Nsim))
        else:
            self.Q = self.InitialList(ones([self.Nsim, self.Nsim], 'd')
                                      - identity(self.Nsim))


    def GetTools(self):
        if self.ens.config.concentrations == "peak":
            return self.Q[0]
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            return self.Q[hour]


    def UpdateTools(self, Q):
        if self.ens.config.concentrations == "peak":
            self.Q[0] = Q
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            self.Q[hour] = Q


    def compute_proba_on_couples(self, Q, loss, prev_weight):
        """
        Computes the matrix of probability on couples (i, j). Its diagonal is
        null.

        @type Q: array
        @type loss: array
        @type prev_weight: array
        """
        N = Q.shape[0]
        # Probability on couples (i, j) with i != j.
        for i in range(N):
            for j in range(N):
                p_ij = dot( i2j(N, i,j), prev_weight)
                lossij = 2. * dot( p_ij, loss )
                Q[i,j] = self.__core_calculus( Q[i,j], lossij)
        Q /= Q.sum()
        return Q


    def __core_calculus(self, prev_w, loss):
        """
        Core function to compute the updated weight vector.

        @type prev_w: array
        @type loss: array
        """
        # Example of EG update.
        return prev_w * exp(- 2. * loss)


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()
        Q = self.GetTools()
        loss = 2. * sum((dot(previous_weight, s) - o) * s, 1)

        Q = self.compute_proba_on_couples(Q, loss, previous_weight)
        weight = self.fixed_point_method(Q, self.precision, previous_weight)
        weight /= weight.sum()

        self.UpdateTools(Q)
        self.AcquireWeight(weight)


#################
# INTERNAL ZINK #
#################


class InternalZink(InternalMethod):
    """
    This class implements the greedy projection gradient algorithm with
    internal regret (Zinkevich, 2003). You should not change the projection
    function ('projection_simplex').
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, extended = False, U = 1., Nskip = 1,
                 Nlearning = 1, learning_rate = 1.e-6,
                 projection_function = projection_simplex, option = "step",
                 precision = 1.e-10, verbose = False, fixed_point_method =
                 fixed_point2):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type learning_rate: float
        @param learning_rate: Learning rate.
        @type projection_function: function
        @param projection_function: Projection on simplex or ball
        """
        self.learning_rate = learning_rate
        self.projection_function = projection_function
        InternalMethod.__init__(self, ens,
                                configuration_file = configuration_file,
                                process = process, statistics = statistics,
                                Nskip = Nskip, Nlearning = Nlearning,
                                option = option, extended = extended,
                                U = U, precision = precision,
                                verbose = verbose,
                                fixed_point_method = fixed_point_method)


    def mat2vec(self, matrix):
        """
        Reshapes a matrix with a null diagonal to a vector (of size n(n-1)).

        @type matrix: 2D array
        @rtype: 1D array
        """
        N = matrix.shape[0]
        out = matrix.reshape((N * N,))
        out = list(out)
        for i in range(N - 1 , -1, -1):
            del out[i * (N + 1)]
        return array(out)


    def vec2mat(self, vec, N):
        """
        Reverses 'mat2vec'. It rebuilds the corresponding matrix with null
        diagonal.

        @type vec: array
        @type N: int
        @rtype: 2D array
        """
        out = list(vec)
        for i in range(N):
            out.insert(i * (N + 1), 0.)
        out = array(out)
        out = out.reshape((N,N))
        return out


    def compute_proba_on_couples(self, Q, loss, prev_weight):
        N = Q.shape[0]
        internal_loss = zeros([N, N], 'd')
        for i in range(N):
            for j in range(N):
                p_ij = dot(i2j(N, i,j), prev_weight)
                lossij = dot(p_ij, loss)
                internal_loss[i, j] = lossij

        Q -= self.learning_rate * internal_loss
        V = self.mat2vec(Q)
        proj_V = self.projection_function(V)
        Q = self.vec2mat(proj_V, N)
        return Q


#######################################
# INTERNAL POLYNOMIAL GRADIENT METHOD #
#######################################


class InternalPolynomialGradient(InternalMethod):
    """
    This class implements polynomially weighted average forecaster using the
    gradient of the loss (Cesa-Bianchi & Lugosi, 2006, pp. 12 & 23) and with
    internal regret.
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, extended = False, U = 1., Nskip = 1,
                 Nlearning = 1, power = 2., option = "step", precision =
                 1.e-10, verbose = False, fixed_point_method = fixed_point2):

        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.
        """
        self.power = power
        InternalMethod.__init__(self, ens,
                                configuration_file = configuration_file,
                                process = process, statistics = statistics,
                                Nskip = Nskip, Nlearning = Nlearning,
                                option = option, extended = extended,
                                U = U, precision = precision,
                                verbose = verbose,
                                fixed_point_method = fixed_point_method)


    def compute_proba_on_couples(self, Q, gradloss, prev_weight):
        """
        Q is the matrix of regret with forecast shifted by index on i and j.

        @type Q: array
        @type gradloss: array
        @type prev_weight: array

        @rtype: array, array
        """
        N = Q.shape[0]
        # regret on couples (i<>j)
        for i in range(N):
            for j in range(N):
                p_ij = prev_weight[i] * (gradloss[i] - gradloss[j])
                Q[i,j] += p_ij
        M = ((Q + abs(Q)) / 2) ** (self.power - 1)
        M /= M.sum()
        return M, Q


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()
        Q = self.GetTools()

        gradloss = 2. * sum((dot(previous_weight, s) - o) * s, 1)
        M, Q = self.compute_proba_on_couples(Q, gradloss, previous_weight)
        weight = self.fixed_point_method(M, self.precision, previous_weight)
        # should not be necessary: postprocessing
        weight /= weight.sum()

        self.UpdateTools(Q)
        self.AcquireWeight(weight)


##############################################
# INTERNAL EXPONENTIATED GRADIENT DISCOUNTED #
##############################################


class InternalExponentiatedGradientDiscounted(InternalMethod):
    """
    This class implements the exponentiated gradient algorithm with discount
    losses (Mallet, Mauricette, and Stoltz, 2007) and with internal regret.
    """

    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, extended = False, U = 1., Nskip = 1,
                 Nlearning = 1, learning_rate = 1.2e-4, forget_rate = 1.,
                 p1 = 0.5, p2 = 1., option = "step", precision = 1.e-10,
                 fixed_point_method = fixed_point2,
                 verbose = False):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments.

        @type learning_rate: float
        @param learning_rate: Learning rate.
        @type forget_rate: float
        @param forget_rate: Forget rate of past confidence.
        @type p1: float
        @param p1: Power of the learning rate.
        @type p2: float
        @param p2: Power of the forget rate.
        """
        self.learning_rate = learning_rate
        self.forget_rate = forget_rate
        self.p1 = p1
        self.p2 = p2
        InternalMethod.__init__(self, ens,
                                configuration_file = configuration_file,
                                process = process, statistics = statistics,
                                Nskip = Nskip, Nlearning = Nlearning,
                                option = option, extended = extended,
                                U = U, precision = precision,
                                verbose = verbose,
                                fixed_point_method = fixed_point_method)


    def Init(self):
        self.initial_weight = ones(self.Nsim, 'd') / float(self.Nsim)
        self.confidence = self.InitialList([])


    def GetTools(self):
        if self.ens.config.concentrations == "peak":
            return self.confidence[0]
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            return self.confidence[hour]


    def UpdateTools(self, confidence):
        if self.ens.config.concentrations == "peak":
            self.confidence[0] = confidence
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            self.confidence[hour] = confidence


    def compute_proba_on_couples(self, confidence):
        """
        Note: this does not work in hourly mode.

        @type confidence: list of array
        @rtype: array
        """
        if self.extended:
            N = 2 * self.Nsim
            if self.ens.config.concentrations == "peak":
                weight_list = self.weight_ext[:]
            elif self.ens.config.concentrations == "hourly":
                raise Exception, "hourly mode not available"
            weight_list.insert(0, self.GetInitialWeight().copy())
        else:
            N = self.Nsim
            if self.ens.config.concentrations == "peak":
                weight_list = self.weight[:]
            elif self.ens.config.concentrations == "hourly":
                raise Exception, "hourly mode not available"
            weight_list.insert(0, self.GetInitialWeight().copy())

        T = len(confidence)
        Q = ones([N, N], 'd') - identity(N)

        for t in range(T):
            for i in range(N):
                for j in range(N):
                    Q[i,j] *= exp( self.learning_rate / (T ** self.p1)
                                  * weight_list[t][i]
                                  * (confidence[t][i] - confidence[t][j])
                                  * (1 + self.forget_rate
                                     / ((T - t) ** self.p2)))
            Q /= Q.sum()
        Q /= Q.sum()
        return Q


    def UpdateWeight(self, s, o):
        previous_weight = self.GetPreviousWeight()
        confidence = self.GetTools()

        loss = 2. * sum((dot(previous_weight, s) - o) * s, 1)
        confidence.append(loss)

        Q = self.compute_proba_on_couples(confidence)
        weight = self.fixed_point_method(Q, self.precision, previous_weight)

        self.UpdateTools(confidence)
        self.AcquireWeight(weight)


#############################
# DYNAMIC LINEAR REGRESSION #
#############################


class DynamicLinearRegression(EnsembleMethod):
    """
    This class implements the Dynamic Linear Regression algorithm
    (West and Harrison 1989).
    """


    def __init__(self, ens, configuration_file = None, process = True,
                 statistics = True, Nskip = 1, Nlearning = 1,
                 verbose = False, delta = 1., s = 1.):
        """
        See documentation of 'EnsembleMethod.__init__' for explanations about
        arguments. Initial vector of forecast should be put in argument

        @type delta: float
        @param delta: rate of uncertainty.
        """
        self.delta = delta
        self.s_tmp = s
        EnsembleMethod.__init__(self, ens,
                                configuration_file = configuration_file,
                                process = process, statistics = statistics,
                                Nskip = Nskip, Nlearning = Nlearning,
                                option = "station", verbose = verbose)


    def Init(self):
        self.initial_weight = zeros(self.Nsim, 'd')
        self.R = self.InitialList(identity(self.Nsim, 'd'))
        self.C = self.InitialList(identity(self.Nsim, 'd'))
        self.W = self.InitialList(zeros([self.Nsim, self.Nsim], 'd'))
        self.e = self.InitialList(zeros(self.Nsim, 'd'))
        self.Q = self.InitialList(1.)
        self.s = self.InitialList(self.s_tmp)
        self.d = self.InitialList(1.)
        self.n = self.InitialList(1)


    def UpdateTools(self, R, Q, C, W, e, s, d, n):
        if self.ens.config.concentrations == "peak":
            self.R[0] = R
            self.Q[0] = Q
            self.C[0] = C
            self.W[0] = W
            self.e[0] = e
            self.s[0] = s
            self.d[0] = d
            self.n[0] = n
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            self.R[hour] = R
            self.Q[hour] = Q
            self.C[hour] = C
            self.W[hour] = W
            self.e[hour] = e
            self.s[hour] = s
            self.d[hour] = d
            self.n[hour] = n


    def GetTools(self):
        if self.ens.config.concentrations == "peak":
            return self.R[0], self.Q[0], self.C[0], self.W[0], self.e[0],\
                   self.s[0], self.d[0], self.n[0]
        elif self.ens.config.concentrations == "hourly":
            hour = self.ens.all_dates[self.step].hour
            return self.R[hour], self.Q[hour], self.C[hour], self.W[hour],\
                   self.e[hour], self.s[hour], self.d[hour], self.n[hour]


    def UpdateWeight(self, sim, obs):
        # R, C, W  variables are variance matrix. Q is a [1, 1]-array.
        prev_weight = self.GetPreviousWeight()
        R, Q, C, W, e, s, d, n = self.GetTools()

        n += 1
        # With known update of W, we would have R = C + W
        # Here, we use a discount factor 'delta'.
        R = C / self.delta
        W = R * (1 - self.delta)
        Q = dot(sim.transpose(), dot(C, sim)) + s
        e = (obs - dot(prev_weight, sim))
        d += s * (e ** 2) / Q[0]
        C = ((d / float(n)) / s) * (R - Q * outer(prev_weight, prev_weight))
        s = (d / float(n))
        A = dot(R, sim) / Q[0]

        weight = prev_weight + e * A.reshape(self.Nsim,)

        self.UpdateTools( R, Q, C, W, e, s, d, n)
        self.AcquireWeight(weight)
