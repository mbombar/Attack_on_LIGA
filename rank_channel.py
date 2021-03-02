#!/usr/bin/env python3

from sage.coding.channel import Channel, format_interval
from sage.rings.integer import Integer
from sage.categories.fields import Fields
from sage.misc.prandom import randint
from sage.modules.free_module import span
from sage.matrix.constructor import matrix
from copy import copy

from linear_rank_metric import from_matrix_representation

class StaticRankErrorChannel(Channel):
    r"""
    Channel which adds an error of static rank to each message it transmits.

    The input space and the output space of this channel are the same.

    INPUT:

    - ``space`` -- the space of both input and output


    - ``rank_error`` -- the rank of the error added to each transmitted message
      It can be either an integer of a tuple. If a tuple is passed as
      argument, the rank of the error will be a random integer between the
      two bounds of the tuple.

    - ``relative_field`` -- The field to which the extension is relative.
      If not given, it will default to the prime_subfield of the ambiant space
      base_field.

    EXAMPLES:

    We construct a StaticRankErrorChannel which adds error of rank 2
    to any transmitted message::

        sage: n_err = 2
        sage: Chan = channels.StaticRankErrorChannel(GF(59)^40, n_err)
        sage: Chan
        Channel creating error of rank 2 over Finite Field of size 59,
        of input and output space Vector space of dimension 40
        over Finite Field of size 59

    We can also pass a tuple for the number of errors::

        sage: n_err = (1, 10)
        sage: Chan = channels.StaticRankErrorChannel(GF(59)^40, n_err)
        sage: Chan
        Channel creating error of rank between 1 and 10 over Finite Field
        of size 59, of input and output space Vector space of dimension 40
        over Finite Field of size 59
    """

    def __init__(self, space, rank_error, relative_field = None):
        r"""
        TESTS:

        If the number of errors exceeds the dimension of the input space,
        it will return an error::

            sage: n_err = 42
            sage: Chan = channels.StaticRankErrorChannel(GF(59)^40, n_err)
            Traceback (most recent call last):
            ...
            ValueError: There might be more errors than the dimension of the input space

        If ``relative_field`` is specified and is not a subfield of the base field,
        it will return an error::

            sage: n_err = 2
            sage: Chan = channels.StaticRankErrorChannel(GF(16)^6, n_err, GF(8))
            Traceback (most recent call last):
            ...
            ValueError: Finite Field in z4 of size 2^4 is not an extension of
            Finite Field in z3 of size 2^3

        If ``relative_field`` is specified and is not a field,
        it will return an error::

            sage: n_err = 2
            sage: Chan = channels.StaticRankErrorChannel(GF(16)^6, n_err, GF(4)^2)
            Traceback (most recent call last):
            ...
            ValueError: relative_field must be a Field and
            Vector space of dimension 2 over Finite Field in z2
            of size 2^2 is not.
        """
        if isinstance(rank_error, (Integer, int)):
            rank_error = (rank_error, rank_error)
        if not isinstance(rank_error, (tuple, list)):
            raise ValueError("rank_error must be a tuple, a list, an Integer or a Python int")
        super(StaticRankErrorChannel, self).__init__(space, space)
        if rank_error[1] > space.dimension():
            raise ValueError("There might be more errors than the dimension of the input space")
        self._rank_error = rank_error
        self._base_field = space.base_field()
        if not relative_field:
            self._relative_field = self._base_field.prime_subfield()
        else:
            if not relative_field in Fields():
                raise ValueError("relative_field must be a Field and %s is not." % relative_field)
            if not relative_field.is_subring(self._base_field):
                raise ValueError("%s is not an extension of %s" % (self._base_field, relative_field))
            self._relative_field = relative_field

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: n_err = 2
            sage: Chan = channels.StaticRankErrorChannel(GF(59)^40, n_err)
            sage: Chan
            Channel creating error of rank 2 over Finite Field of size 59,
            of input and output space Vector space of dimension 40
            over Finite Field of size 59
        """
        no_err = self.rank_error()
        return "Channel creating error of rank %s over %s, of input and output space %s"\
                    % (format_interval(no_err), self._relative_field, self.input_space())

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: n_err = 2
            sage: Chan = channels.StaticRankErrorChannel(GF(59)^40, n_err)
            sage: latex(Chan)
            \textnormal{Channel creating error of rank 2 over Finite Field
            of size 59, of input and output space Vector space of dimension 40
            over Finite Field of size 59}
        """
        no_err = self.rank_error()
        return "\\textnormal{Channel creating error of rank %s over %s, of input and output space %s}"\
                % (format_interval(no_err), self._relative_field, self.input_space())

    def transmit_unsafe(self, message):
        r"""
        Returns ``message`` with as many errors as ``self._rank_error`` in it.

        If ``self._rank_error`` was passed as a tuple for the number of errors, it will
        pick a random integer between the bounds of the tuple and use it as the number of errors.

        This method does not check if ``message`` belongs to the input space of``self``.

        INPUT:

        - ``message`` -- a vector

        OUTPUT:

        - a vector of the output space

        EXAMPLES::

            sage: F = GF(16)^6
            sage: n_err = 2
            sage: Chan = channels.StaticRankErrorChannel(F, n_err, GF(4))
            sage: set_random_seed(10)
            sage: msg = F.random_element()
            sage: msg
            (z4 + 1, z4, z4^3 + z4 + 1, z4^3 + z4^2 + z4 + 1, z4^2, z4^2)
            sage: set_random_seed(10)
            sage: c = Chan.transmit_unsafe(msg)
            (z4^3 + z4, 1, z4^2 + 1, z4^3 + z4^2 + z4 + 1, 1, z4^3 + z4^2 + 1)

        TESTS::

            sage: from sage.coding.linear_rank_metric import rank_distance
            sage: rank_distance(msg, c, GF(4))
            2
        """
        w = copy(message)
        rank_error = randint(*self.rank_error())
        Fqm = self._base_field
        Fq = self._relative_field
        V = Fqm.vector_space(Fq, map=False)
        n = self.input_space().dimension()
        good = False
        w = None
        while not good:
            basis = [V.random_element() for i in range(rank_error)]
            R = span(basis)
            err = [R.random_element() for i in range(n)]
            M = matrix(Fq, err).transpose()
            if M.rank() == rank_error:
                good = True
        e = from_matrix_representation(M, Fqm)
        w = message + e
        return w

    def rank_error(self):
        r"""
        Returns the number of errors created by ``self``.

        EXAMPLES::

            sage: n_err = 3
            sage: Chan = channels.StaticRankErrorChannel(GF(59)^6, n_err)
            sage: Chan.rank_errors()
            (3, 3)
        """
        return self._rank_error

    def number_errors(self):
        r"""
        This function is here to ease the life of people comming from the Hamming world.
        Returns the number of errors created by ``self``.

        EXAMPLES::

            sage: n_err = 3
            sage: Chan = channels.StaticRankErrorChannel(GF(59)^6, n_err)
            sage: Chan.number_errors()
            (3, 3)
        """
        return self._rank_error
