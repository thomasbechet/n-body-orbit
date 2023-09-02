import numpy as np
import matplotlib.pyplot as plt

G = 6.67408e-11 # m^3 kg^-1 s^-2

class Orbit:
    def __init__(self, a, e, om, i, w, T):
        self.a = a
        self.e = e
        self.om = om
        self.i = i
        self.w = w
        self.T = T

    def cartesian(self, t0, t):

        # T = 2.0 * np.pi * np.sqrt((self.a**3) / mu) # Orbital period
        if t == t0:
            t = t0
            Mt = 0
        else:
            deltaT = 2.0 * np.pi * ((t - t0) / self.T)
            Mt = deltaT

        # Newton's method
        E = Mt
        F = E - self.e * np.sin(E) - Mt
        j = 0
        while np.abs(F) > 0.00001 and j < 30:
            E = E - F / (1.0 - self.e * np.cos(E))
            F = E - self.e * np.sin(E) - Mt
            j += 1

        # True anomaly
        nu = 2.0 * np.arctan2(np.sqrt(1.0 + self.e) * np.sin(E / 2.0), np.sqrt(1.0 - self.e) * np.cos(E / 2.0))

        # Radius to central body
        rc = self.a *(1.0 - self.e * np.cos(E))

        X = rc * (np.cos(self.om) * np.cos(self.w + nu) - np.sin(self.om) * np.sin(self.w + nu) * np.cos(self.i))
        Y = rc * (np.sin(self.om) * np.cos(self.w + nu) + np.cos(self.om) * np.sin(self.w + nu) * np.cos(self.i))
        Z = rc * (np.sin(self.i) * np.sin(self.w + nu))

        return np.array([X, Y, Z])

    def draw(self, ax, refpos):
        times = np.linspace(0.0, self.T, 360)
        x = []
        y = []
        z = []
        for t in times:
            pos = self.cartesian(0.0, t)
            x.append(refpos[0] + pos[0])
            y.append(refpos[1] + pos[1])
            z.append(refpos[2] + pos[2])
        x.append(refpos[0])
        y.append(refpos[1])
        z.append(refpos[2])
        ax.plot(x, y, z, zorder=10, c='b')

class Planet:

    def __init__(self, mass, radius, pos):
        self.mass = mass
        self.radius = radius
        self.pos = pos
        self.mu = G * mass

    def compute_force(self, pos):
        rel = self.pos - pos
        d = np.linalg.norm(rel)
        return (self.mu / (d * d)) * (rel / d)

    def draw(self, ax):
        u, v = np.mgrid[0:2 * np.pi:30j, 0:np.pi:20j]
        x = self.pos[0] + self.radius * np.cos(u) * np.sin(v)
        y = self.pos[1] + self.radius * np.sin(u) * np.sin(v)
        z = self.pos[2] + self.radius * np.cos(v)
        ax.plot_surface(x, y, z, zorder=1, linewidth=0, antialiased=False)

class Body:

    def __init__(self, pos, vel):
        self.pos = pos
        self.vel = vel
        self.acc = np.array([0.0, 0.0, 0.0], dtype=np.float32)

    def update(self, planets, dt):
        # New position
        npos = (dt * dt * 0.5 * self.acc) + (self.vel * dt) + self.pos
        # New acceleration
        nacc = 0.0
        for planet in planets.values():
            nacc += planet.compute_force(npos)
        # New velocity
        nvel = (nacc + self.acc) * dt * 0.5 + self.vel
        # Update body
        self.pos = npos
        self.vel = nvel
        self.acc = nacc

    def energy(self):
        v = 0.5 * np.linalg.norm(self.vel)**2
        rel = epos - self.pos
        d = np.linalg.norm(rel)
        p = -((G * Me) / d)
        return v + p

def days_to_seconds(days):
    return days * (60.0 * 60.0 * 24.0)
def seconds_to_days(secs):
    return secs / (60.0 * 60.0 * 24.0)

def simulate(simu):
    import json
    # Load configuration
    f = open('sim.json')
    config = json.load(f)
    f.close()
    # Get system and simulation
    sim = config['simulations'][simu]
    system = config['systems'][sim['system']]
    # Load system
    planets = {}
    orbits = {}
    for name, info in system.items():
        if name.startswith('_'):
            continue
        if info['type'] == 'fixed':
            planet = Planet(info['mass'], info['radius'], np.array(info['pos']))
            planets[name] = planet
        elif info['type'] == 'orbit':
            planet = Planet(info['mass'], info['radius'], np.array([0.0, 0.0, 0.0]))
            planets[name] = planet
            oinfo = info['orbit']
            a = 0.5
            e = 0.0
            if 'apo' in oinfo and 'per' in oinfo:
                apo = oinfo['apo']
                per = oinfo['per']
                a = (apo + per) / 2
                e = (apo - per) / (apo + per)
            else:
                a = oinfo['a']
                e = oinfo['e']
            om = np.radians(oinfo['om'])
            i = np.radians(oinfo['i'])
            w = np.radians(oinfo['w'])
            T = oinfo['T']
            orbit = Orbit(a, e, om, i, w, days_to_seconds(T))
            orbits[name] = {
                    'ref': oinfo['ref'],
                    'orbit': orbit,
                    }

    # Load bodies
    bodies = {}
    for name, info in sim['bodies'].items():
        body = Body(np.array(info['pos']), np.array(info['vel']))
        bodies[name] = body

    # Load simulation info
    simtime = days_to_seconds(sim['simtime'])
    simdelta = sim['simdelta']
    simstep = int(simtime / simdelta)
    print('Days:', seconds_to_days(simtime))
    print('TimeStep:', simdelta)
    print('Steps:', simstep)

    #Prepare recording
    bodies_rec = {}
    for name, body in bodies.items():
        relatives = {}
        bodies_rec[name] = {
                'x_pos': [],
                'y_pos': [],
                'z_pos': [],
                'crash': False,
                'relatives': {},
                'color': sim['bodies'][name]['color']
                }

    ### START SIMULATION ###
    i = 0
    t = 0.0
    while i < simstep:
        # Update orbit planets
        for name, orbit in orbits.items():
            planet = planets[name]
            kep = orbit['orbit']
            ref_planet = planets[orbit['ref']]
            planet.pos = ref_planet.pos + kep.cartesian(days_to_seconds(-9.2), t)

        # Update bodies
        for name, body in bodies.items():
            # Check crashed
            if bodies_rec[name]['crash']:
                continue

            # Update
            body.update(planets, simdelta)

            # Detect crash
            for _, planet in planets.items():
                if np.linalg.norm(planet.pos - body.pos) <= planet.radius:
                    bodies_rec[name]['crash'] = True
                    print(name, 'crashed on day:', seconds_to_days(t))
                    break
            if not bodies_rec[name]['crash']:
                if i % 4 == 0:
                    bodies_rec[name]['x_pos'].append(body.pos[0])
                    bodies_rec[name]['y_pos'].append(body.pos[1])
                    bodies_rec[name]['z_pos'].append(body.pos[2])

            # Print progression
            if i % 10000 == 0:
                print('{:.2f}%'.format((t / simtime) * 100.0))

            # Next step
            i += 1
            t += simdelta
    ### END SIMULATION ###

    ### DISPLAY FINAL SYSTEM STATE ###
    from matplotlib import cm
    fig = plt.figure()
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect(aspect = (1, 1, 1))
    ax.set_xlim((-1e8, 1e8))
    ax.set_ylim((-1e8, 1e8))
    ax.set_zlim((-1e8, 1e8))
    for planet in planets.values():
        planet.draw(ax)
    for name, orbit in orbits.items():
        pos = planets[name].pos
        ref = planets[orbit['ref']].pos
        orbit['orbit'].draw(ax, ref)
        ax.plot(
                (pos[0], ref[0]),
                (pos[1], ref[1]),
                (pos[2], ref[2]),
                zorder=10, c='g'
                )
    for rec in bodies_rec.values():
        color = rec['color']
        plt.plot(rec['x_pos'], rec['y_pos'], rec['z_pos'], c=color, zorder=10)

    # Draw relative trajectory
    #rel_x_pos = [planets[1].pos[0] + rel[0] for rel in rel_moon_pos]
    #rel_y_pos = [planets[1].pos[1] + rel[1] for rel in rel_moon_pos]
    #rel_z_pos = [planets[1].pos[2] + rel[2] for rel in rel_moon_pos]
    #plt.plot(rel_x_pos, rel_y_pos, rel_z_pos, c='m', zorder=10)
    #times = np.linspace(0.0, nstep * sim_dt, num=nstep)

    plt.show()

import argparse

parser = argparse.ArgumentParser(description='N-Body simulation')
parser.add_argument('simulation', metavar='simulation', type=str, help='Simulation to run')
args = parser.parse_args()
simulate(args.simulation)
