from sage.all import *
from quiver import *

class QuiverSetting():
    
    def __init__(self, Q, d) -> None:
        self._Q = Q
        self._d = d

    def sst_inequalities(self):
        r"""
        Returns the inequalities for the sst cone of `d`.
        """

        r"""The inequalities are given as a list of lists; each such list [`a_1`,...,`a_n`] is to be understood as the inequality :math:`a_1x_1 + ... + a_nx_n \geq 0`."""

        Q, d = self._Q, self._d
        genSubdims = Q.all_generic_subdimension_vectors(d, proper=True, nonzero=True)
        ineqs = [list(-e) for e in genSubdims]
        return ineqs
    
    def sst_cone(self):
        r"""
        Returns the sst cone of `d`.
        """

        r"""It is given by the equality :math:`\theta(d) = 0` and the inequalities from :meth:`sst_inequalities()`."""

        d = self._d 
        ineqs = self.sst_inequalities()
        ieqs = [[0]+ineq for ineq in ineqs]
        return Polyhedron(ieqs=ieqs, eqns=[[0]+list(d)])
    
    def hyperplane_cap_sst_cone(self, e):
        r"""
        Returns the intersection of hyperplane math:`H_e` with math:`\operatorname{sst}(d)`.
        """

        d = self._d
        ineqs = self.sst_inequalities()
        ieqs = [[0]+ineq for ineq in ineqs]
        return Polyhedron(ieqs=ieqs, eqns=[[0]+list(d), [0]+list(e)])
    
    def wall_inequalities(self, e):
        r"""
        Returns the inequalities for the wall :math:`W_e`.
        """

        r""":math:`W_e` is the intersection :math:`\mathrm{sst}(e) \cap \mathrm{sst}(d-e)`."""

        Q, d = self._Q, self._d
        Qe = QuiverSetting(Q, e)
        Qf = QuiverSetting(Q, d-e)
        return Qe.sst_inequalities() + Qf.sst_inequalities()
    
    def wall(self, e):
        r"""
        Returns the wall :math:`W_e`.
        """

        d = self._d
        ineqs = self.wall_inequalities(e)
        ineqs = [[0]+ineq for ineq in ineqs]
        return Polyhedron(ieqs=ineqs, eqns=[[0]+list(e), [0]+list(d-e)])
    
    def wall_contains_hyperplane(self, e):
        r"""
        Returns true if the intersection of :math:`H_e` with the semistable cone is contained in the wall :math:`W_e`
        """

        W = self.wall(e)
        H = self.hyperplane_cap_sst_cone(e)
        return all([W.contains(list(r)) for r in H.ray_generator()])
    
    def all_walls(self):
        r"""
        Returns the list of dicts of all walls.
        """

        Q, d = self._Q, self._d
        subdims = Q.all_subdimension_vectors(d, proper=True, nonzero=True)
        walls = []
        for e in subdims:
            W = self.wall(e)
            contained = False
            for j in range(len(walls)):
                Wdatum = walls[j]
                if W == Wdatum["wall"]:
                    contained = True
                    Wdatum["subdimension"] += [e]
                    break
            if not contained:
                Wdatum = {
                    "wall" : self.wall(e),
                    "subdimension" : [e],
                    "inequalities" : self.wall_inequalities(e)
                }
                walls += [Wdatum]

        return walls
    
    def git_equivalent(self, theta, eta, walls=None):
        r"""
        Returns the truth value of `theta` and `eta` being GIT equivalent.
        """

        r"""Two stability parameters :math:`\theta` and :math:`\eta` are GIT equivalent iff for every wall :math:`W_e` holds that the convex hull :math:`[\theta, \eta]` is contained in :math:`W_e` or does not intersect it."""

        equivalent = True 
        Q, d = self._Q, self._d
        if walls == None:
            walls = self.all_walls()
        subdims = [W["subdimension"][0] for W in walls]
        ineqsList = [W["inequalities"] for W in walls]
        for i in range(len(subdims)):
            e = subdims[i]
            ineqs = ineqsList[i]
            # Check if the line segment [theta, eta] is parallel to the hyperplane spanned by W_e 
            parallel = ((eta-theta)*e == 0) 
            if parallel:
                # Check for containment of [theta, eta] in W_e
                # As W_e is convex, it's enough to check if theta and eta are contained
                equivalent = equivalent and all([vector(ineq)*theta >= 0 and vector(ineq)*eta >= 0 for ineq in ineqs])
            else: 
                # In this case the line through theta and eta and the hyperplane H_e spanned by W_e are transversal
                # Compute the unique intersection of the line through theta and eta and the hyperplane H_e 
                t = -eta*e/((theta-eta)*e)
                # If t notin [0,1] then [theta, eta] doesn't intersect H_e (and thus also not W_e)
                disjoint = (t < 0 or t > 1)
                if not disjoint:
                    zeta = t*theta + (1-t)*eta
                    # If zeta doesn't satisfy all the inequalities of W_e, then [theta, eta] is disjoint from W_e
                    disjoint = any([vector(ineq)*zeta < 0 for ineq in ineqs])
                equivalent = equivalent and disjoint
            if not equivalent:
                break

        return equivalent
    
    def is_generic(self, theta, walls=None):
        if walls == None:
            walls = self.all_walls()
        for Wdatum in walls:
            W = Wdatum["wall"]
            if W.contains(theta):
                return False
        return True
    
    def render_2d(self, projectionAxis=None):
        r"""
        Plots the walls for a quiver with 3 vertices.
        """

        Q, d = self._Q, self._d 
        assert Q.number_of_vertices() == 3
        if projectionAxis == None:
            i = 0 # Project along first axis by default
        else:
            i = projectionAxis
        assert i in Q.vertices() and d[i] != 0

        walls3d = [W["wall"] for W in self.all_walls()] # Just the polyhedra
        verticesList = [W.vertices_list() for W in walls3d]
        linesList = [W.lines_list() for W in walls3d]
        raysList = [W.rays_list() for W in walls3d]

        # Project down to two dimensions
        verticesList = [[vector([vertex[j] for j in range(3) if j != i]) for vertex in vertices] for vertices in verticesList]        
        linesList = [[vector([line[j] for j in range(3) if j != i]) for line in lines] for lines in linesList]
        raysList = [[vector([ray[j] for j in range(3) if j != i]) for ray in rays] for rays in raysList]

        walls2d = [Polyhedron(vertices=verticesList[j], lines=linesList[j], rays=raysList[j]) for j in range(len(verticesList))]

        return sum([walls2d[j].plot() for j in range(len(walls2d))])
    
    def random_stability_parameter(self, min=-1, max=1, requireSst=False, requireGeneric=False, walls=None):
        r"""
        Generates a random stability parameter for :math:`Q`, :math:`d`.

        If requireSst is `True` then it continues searching until :math:`d` is semistable for :math:`theta`.
        If requireGeneric is `True` then it continues searching until :math:`theta` does not lie on a wall.
        """

        # TODO: Error handling not yet implemented.

        Q = self._Q
        d = self._d
        m = Q.number_of_vertices()

        assert all([d[i] != 0 for i in range(m)])

        while True:
            x = _random_rational_vector(m-1, min=min, max=max)
            y = -sum([x[i]*d[i] for i in range(m-1)])/d[m-1]
            theta = vector(x+[y])

            if not requireSst or Q.has_semistable_representation(d, theta):
                break
            if not requireGeneric or self.is_generic(theta, walls=walls):
                break

        return theta

"""
Helper functions
"""

def _random_rational_vector(m, min=-1, max=1, error=0.001):
    return [n(uniform(min, max)).nearby_rational(max_error=error) for _ in range(m)]
