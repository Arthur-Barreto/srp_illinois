import numpy as np

def eebls(n, t, x, nf, fmin, df, nb, qmi, qma):
    """
    Compute the BLS spectrum.

    Parameters:
    n    : int     - Number of data points
    t    : array   - Array of time values
    x    : array   - Array of data values
    nf   : int     - Number of frequency points in the spectrum
    fmin : float   - Minimum frequency (must be > 0)
    df   : float   - Frequency step
    nb   : int     - Number of bins in the folded time series
    qmi  : float   - Minimum fractional transit length
    qma  : float   - Maximum fractional transit length

    Returns:
    p    : array   - Values of the BLS spectrum
    bper : float   - Period at the highest peak in the spectrum
    bpow : float   - Value of the BLS spectrum at the highest peak
    depth: float   - Depth of the transit at *bper*
    qtran: float   - Fractional transit length
    in1  : int     - Bin index at the start of the transit
    in2  : int     - Bin index at the end of the transit
    """

    minbin = 5
    nbmax = 2000

    if nb > nbmax:
        raise ValueError("nb > nbmax")

    tot = t[n-1] - t[0]
    if fmin < 1.0 / tot:
        raise ValueError("fmin < 1/total time span")

    rn = float(n)
    kmi = int(qmi * nb)
    if kmi < 1:
        kmi = 1

    kma = int(qma * nb) + 1
    kkmi = int(rn * qmi)
    if kkmi < minbin:
        kkmi = minbin

    bpow = 0.0

    u = t - t[0]
    s = np.mean(x)
    v = x - s

    p = np.zeros(nf)

    for jf in range(nf):
        f0 = fmin + df * jf
        p0 = 1.0 / f0

        y = np.zeros(nb)
        ibi = np.zeros(nb, dtype=int)

        for i in range(n):
            ph = u[i] * f0
            ph -= int(ph)
            j = int(nb * ph)
            ibi[j] += 1
            y[j] += v[i]

        y = np.concatenate((y, y[:kma]))
        ibi = np.concatenate((ibi, ibi[:kma]))

        power = 0.0

        for i in range(nb):
            s = 0.0
            k = 0
            kk = 0
            nb2 = i + kma

            for j in range(i, nb2):
                k += 1
                kk += ibi[j]
                s += y[j]
                if k < kmi or kk < kkmi:
                    continue
                rn1 = float(kk)
                pow = s * s / (rn1 * (rn - rn1))
                if pow > power:
                    power = pow
                    jn1 = i
                    jn2 = j
                    rn3 = rn1
                    s3 = s

        power = np.sqrt(power)
        p[jf] = power

        if power > bpow:
            bpow = power
            in1 = jn1
            in2 = jn2
            qtran = rn3 / rn
            depth = -s3 * rn / (rn3 * (rn - rn3))
            bper = p0

    if in2 > nb:
        in2 -= nb

    return p, bper, bpow, depth, qtran, in1, in2

# Example usage
n = 100  # number of data points
t = np.linspace(0, 100, n)  # time array
x = np.random.randn(n)  # data array with random noise
nf = 1000  # number of frequency points
fmin = 0.01  # minimum frequency
df = 0.001  # frequency step
nb = 100  # number of bins
qmi = 0.01  # minimum fractional transit length
qma = 0.1  # maximum fractional transit length

p, bper, bpow, depth, qtran, in1, in2 = eebls(n, t, x, nf, fmin, df, nb, qmi, qma)

print("BLS Spectrum:", p)
print("Best Period:", bper)
print("Best Power:", bpow)
print("Depth:", depth)
print("Fractional Transit Length:", qtran)
print("Transit Start Bin:", in1)
print("Transit End Bin:", in2)
