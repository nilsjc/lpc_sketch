namespace LPCcoder
{
    using NAudio.Wave;
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Numerics;
    using System.Globalization;
    internal class Program
    {
        static void Main(string[] args)
        {
            LPC lpc = new();
            Console.WriteLine("Enter the path to your wav file:");
            string filePath = Console.ReadLine()?.Trim('"');
            if(filePath == null || !System.IO.File.Exists(filePath))
            {
                filePath = "voice.wav";
                Console.WriteLine($"File not found. Using default: {filePath}");
            }
            Console.WriteLine("Robot voice?");
            bool useFixedPitch = Console.ReadLine()?.Trim().ToLower() == "y";
            if(useFixedPitch)
            {
                Console.WriteLine("Please enter frequency of the pitch in Hz (e.g. 70)");
                if (int.TryParse(Console.ReadLine(), out int pitchHz))
                    lpc.FixedPitchHz = pitchHz;
            }
            Console.WriteLine("Formant multiplier (1.0 = unchanged, > 1.0 = brighter/smaller voice, < 1.0 = darker/larger, e.g. 1.15 or 0.85):");
            if (float.TryParse(Console.ReadLine(), NumberStyles.Float, CultureInfo.InvariantCulture, out float formantScale))
                lpc.FormantScale = formantScale;
            else
                Console.WriteLine("Invalid input for formant scale. Using default value of 1.0 (no change).");

            Console.WriteLine("Pitch modulation in span like -300 to 300 (0 = none, positive = higher, negative = lower):");
            if (int.TryParse(Console.ReadLine(), out int pitchMod))
                lpc.PitchModulation = pitchMod;

            lpc.PerformLPCAnalysisSynthesizing(filePath, useFixedPitch: useFixedPitch);
            Console.WriteLine("Press any key to exit.");
            Console.ReadKey();
        }
    }
 
    /// <summary>Per-frame analysis result.</summary>
    public class FrameInfo
    {
        public float[] Lpc = Array.Empty<float>(); // a_1 .. a_Order (A(z) = 1 + sum a_k z^-k)
        public float Gain;                          // per-sample excitation amplitude
        public bool Voiced;
        public int PitchPeriod;                     // in samples (valid when Voiced)
    }
 
    public class LPC
    {
        private const int Order = 20;
        private const float PreEmphasis = 0.97f;
        private const float FrameSeconds = 0.030f;  // 30 ms analysis window
        private const float VoicingThreshold = 0.30f;
 
        // TODO not functioning well
        // Formant shift: multiplies every formant frequency by this factor.
        // 1.0 = unchanged, > 1.0 = brighter/smaller voice, < 1.0 = darker/larger.
        // Try 1.15 or 0.85 to hear the effect.
        public float FormantScale = 1.0f;

        // When using a constant pitch instead of the estimated one
        public int FixedPitchHz = 70;
        public int PitchModulation = 0;
 
 
        public void PerformLPCAnalysisSynthesizing(string filePath, bool useFixedPitch = false)
        {
            int sampleRate;
            float[] rawSignal = ReadWav(filePath, out sampleRate);
            float[] signal = CleanUp(rawSignal);
            signal = PreEmphasize(signal, PreEmphasis);
 
            int frameLength = (int)(FrameSeconds * sampleRate);
            int hop = frameLength / 2;               // 50 % overlap for analysis
 
            List<float[]> frames = CreateFrames(signal, frameLength, hop);
 
            var frameData = new List<FrameInfo>(frames.Count);
            foreach (var fr in frames)
                frameData.Add(AnalyzeFrame(fr, Order, sampleRate));
 
            Synthesize(frameData, hop, sampleRate, "Result.wav", useFixedPitch);
            Console.WriteLine("Your soundfile is ready");
        }
 
        /// <summary>
        /// Reads a wav file into 32-bit float samples (mono mix as provided by NAudio).
        /// </summary>
        public float[] ReadWav(string filename, out int sampleRate)
        {
            using var reader = new AudioFileReader(filename);
            sampleRate = reader.WaveFormat.SampleRate;
 
            // AudioFileReader exposes IEEE float samples, so sample count = bytes / 4.
            int bytesPerSample = reader.WaveFormat.BitsPerSample / 8; // = 4 for float
            int capacity = (int)(reader.Length / bytesPerSample);
 
            var buffer = new float[capacity];
            int total = 0;
            int read;
            // Read until the stream is exhausted; honour the returned count.
            while (total < buffer.Length &&
                   (read = reader.Read(buffer, total, buffer.Length - total)) > 0)
            {
                total += read;
            }
 
            if (total != buffer.Length)
                Array.Resize(ref buffer, total);
 
            // If the source is stereo, down-mix to mono.
            int channels = reader.WaveFormat.Channels;
            if (channels > 1)
            {
                int frames = buffer.Length / channels;
                var mono = new float[frames];
                for (int i = 0; i < frames; i++)
                {
                    float s = 0f;
                    for (int c = 0; c < channels; c++) s += buffer[i * channels + c];
                    mono[i] = s / channels;
                }
                buffer = mono;
            }
 
            return buffer;
        }
 
        /// <summary>Trims leading silence (zero samples) from the signal.</summary>
        public float[] CleanUp(float[] input)
        {
            int index = 0;
            while (index < input.Length && input[index] == 0.0f) index++;
            return input.Skip(index).ToArray();
        }
 
        /// <summary>Pre-emphasis filter y[n] = x[n] - a*x[n-1] (boosts high formants).</summary>
        public float[] PreEmphasize(float[] input, float coeff)
        {
            if (input.Length == 0) return input;
            var output = new float[input.Length];
            output[0] = input[0];
            for (int n = 1; n < input.Length; n++)
                output[n] = input[n] - coeff * input[n - 1];
            return output;
        }
 
        /// <summary>
        /// Splits the signal into overlapping, Hann-windowed frames.
        /// </summary>
        public List<float[]> CreateFrames(float[] input, int frameSize, int hop)
        {
            var result = new List<float[]>();
            if (input.Length < frameSize) return result;
 
            for (int start = 0; start + frameSize <= input.Length; start += hop)
            {
                var frame = new float[frameSize];
                for (int x = 0; x < frameSize; x++)
                    frame[x] = input[start + x] * CalculateHanning(x, frameSize);
                result.Add(frame);
            }
            return result;
        }
 
        /// <summary>Hann window value at a given position (reduces spectral leakage).</summary>
        public float CalculateHanning(int index, int length)
        {
            return (float)(0.5 * (1 - Math.Cos((2 * Math.PI * index) / (length - 1))));
        }
 
        /// <summary>
        /// Full analysis of one frame: LPC coefficients, excitation gain,
        /// and a voiced/unvoiced decision with pitch estimate.
        /// </summary>
        public FrameInfo AnalyzeFrame(float[] frame, int order, int sampleRate)
        {
            var lpcCoeffs = new float[order + 1];
            var autocorr = new float[order + 1];
            var reflectionCoeffs = new float[order + 1];
 
            float error = PerformLPCAnalysis(frame, order, lpcCoeffs, autocorr, reflectionCoeffs);
 
            var info = new FrameInfo
            {
                // a_1 .. a_order (drop the leading 1.0)
                Lpc = lpcCoeffs.Skip(1).Take(order).ToArray(),
                // residual energy -> per-sample amplitude
                Gain = (float)Math.Sqrt(Math.Max(error, 0f) / frame.Length)
            };

            // Optionally move the formants by scaling the pole angles.
            // Moving the poles also changes the all-pole filter's overall gain
            // (often by several orders of magnitude), which would swamp the timbre
            // change with a huge volume swing. So we measure the filter energy
            // before and after and rescale the excitation gain to keep loudness
            // constant -- only the formants move, not the volume.
            if (Math.Abs(FormantScale - 1.0f) > 1e-6f)
            {
                float energyBefore = ImpulseResponseEnergy(info.Lpc, 2048);
                info.Lpc = ApplyFormantScale(info.Lpc, FormantScale);
                float energyAfter = ImpulseResponseEnergy(info.Lpc, 2048);
                if (energyAfter > 1e-9f)
                    info.Gain *= (float)Math.Sqrt(energyBefore / energyAfter);
            }
 
            EstimatePitch(frame, sampleRate, out bool voiced, out int period);
            info.Voiced = voiced;
            info.PitchPeriod = period;
            return info;
        }
 
        /// <summary>
        /// Autocorrelation + Levinson-Durbin recursion.
        /// Returns the final prediction-error energy.
        /// </summary>
        public static float PerformLPCAnalysis(float[] frame, int order,
            float[] lpcCoeffs, float[] autocorr, float[] reflectionCoeffs)
        {
            int frameSize = frame.Length;
 
            // Step 1: autocorrelation coefficients R[0..order]
            for (int k = 0; k <= order; k++)
            {
                float sum = 0.0f;
                for (int n = 0; n < frameSize - k; n++)
                    sum += frame[n] * frame[n + k];
                autocorr[k] = sum;
            }
 
            // Step 2: Levinson-Durbin recursion
            lpcCoeffs[0] = 1.0f;
            float error = autocorr[0];
 
            if (error <= 0.0f)
            {
                for (int i = 0; i <= order; i++) lpcCoeffs[i] = 0.0f;
                return 0.0f;
            }
 
            reflectionCoeffs[0] = -autocorr[1] / autocorr[0];
            lpcCoeffs[1] = reflectionCoeffs[0];
            error *= (1.0f - reflectionCoeffs[0] * reflectionCoeffs[0]);
 
            for (int m = 1; m < order; m++)
            {
                float sum = 0.0f;
                for (int j = 0; j <= m; j++)
                    sum += lpcCoeffs[j] * autocorr[m + 1 - j];
 
                reflectionCoeffs[m] = (error != 0.0f) ? -sum / error : 0.0f;
 
                for (int j = 1; j <= (m + 1) / 2; j++)
                {
                    float tmp = lpcCoeffs[j] + reflectionCoeffs[m] * lpcCoeffs[m + 1 - j];
                    lpcCoeffs[m + 1 - j] += reflectionCoeffs[m] * lpcCoeffs[j];
                    lpcCoeffs[j] = tmp;
                }
                lpcCoeffs[m + 1] = reflectionCoeffs[m];
                error *= (1.0f - reflectionCoeffs[m] * reflectionCoeffs[m]);
            }
 
            return error;
        }
 
        /// <summary>
        /// Simple pitch estimate via normalized autocorrelation peak in the
        /// 70-350 Hz range. Decides voiced/unvoiced from the peak strength.
        /// </summary>
        private void EstimatePitch(float[] frame, int sampleRate, out bool voiced, out int period)
        {
            int minLag = Math.Max(1, sampleRate / 350);
            int maxLag = Math.Min(frame.Length - 1, sampleRate / 70);
 
            float r0 = 0f;
            for (int n = 0; n < frame.Length; n++) r0 += frame[n] * frame[n];
 
            voiced = false;
            period = sampleRate / 120; // default fallback pitch (~120 Hz)
 
            if (r0 <= 0f) return;
 
            float bestValue = 0f;
            int bestLag = -1;
            for (int lag = minLag; lag <= maxLag; lag++)
            {
                float sum = 0f;
                for (int n = 0; n < frame.Length - lag; n++)
                    sum += frame[n] * frame[n + lag];
 
                float norm = sum / r0;
                if (norm > bestValue)
                {
                    bestValue = norm;
                    bestLag = lag;
                }
            }
 
            if (bestValue >= VoicingThreshold && bestLag > 0)
            {
                voiced = true;
                period = bestLag;
            }
        }
 
        /// <summary>
        /// Shifts every formant by scaling the pole angles by 'alpha'.
        /// The formants are the poles of 1/A(z); a pole at angle theta sits at
        /// frequency f = theta * fs / (2*pi). Multiplying every angle by alpha
        /// therefore multiplies every formant frequency by alpha, while the pole
        /// radius (and thus the formant bandwidth) is preserved.
        /// </summary>
        /// <param name="aCoeffs">LPC coefficients a_1..a_order.</param>
        /// <param name="alpha">Formant scale factor (1.0 = unchanged).</param>
        private float[] ApplyFormantScale(float[] aCoeffs, float alpha)
        {
            int p = aCoeffs.Length;
 
            // Polynomial z^p + a_1 z^(p-1) + ... + a_p, coefficients high->low.
            var poly = new Complex[p + 1];
            poly[0] = Complex.One;
            for (int i = 0; i < p; i++) poly[i + 1] = new Complex(aCoeffs[i], 0.0);
 
            Complex[] roots = FindRoots(poly);
            // Scale each pole's angle; keep its radius (clamp to stay stable).
            const double maxAngle = Math.PI * 0.999;
            for (int i = 0; i < roots.Length; i++)
            {
                double mag = roots[i].Magnitude;
                double ang = roots[i].Phase * alpha;
                if (ang > maxAngle) ang = maxAngle;
                if (ang < -maxAngle) ang = -maxAngle;
                if (mag >= 0.999) mag = 0.999;       // guarantee stability
                roots[i] = Complex.FromPolarCoordinates(mag, ang);
            }
 
            // Rebuild the monic polynomial: product of (z - root_i).
            var newPoly = new Complex[] { Complex.One };
            foreach (var r in roots)
            {
                var next = new Complex[newPoly.Length + 1];
                for (int i = 0; i < newPoly.Length; i++)
                {
                    next[i]     += newPoly[i];        // multiply by z
                    next[i + 1] += newPoly[i] * (-r); // multiply by (-root)
                }
                newPoly = next;
            }
 
            // Imaginary parts cancel for conjugate pole pairs; keep the real part.
            var result = new float[p];
            for (int i = 0; i < p; i++) result[i] = (float)newPoly[i + 1].Real;
            return result;
        }
 
        /// <summary>
        /// Energy of the all-pole filter's impulse response, i.e. sum of h[n]^2
        /// for H(z) = 1 / A(z), A(z) = 1 + sum aCoeffs[k] z^-(k+1).
        /// Used to keep loudness constant when the formants are moved.
        /// </summary>
        private float ImpulseResponseEnergy(float[] aCoeffs, int n)
        {
            int p = aCoeffs.Length;
            var h = new float[n];
            float energy = 0f;
            for (int i = 0; i < n; i++)
            {
                float acc = (i == 0) ? 1.0f : 0.0f;   // unit impulse input
                for (int k = 0; k < p; k++)
                    if (i - 1 - k >= 0) acc -= aCoeffs[k] * h[i - 1 - k];
                h[i] = acc;
                energy += acc * acc;
            }
            return energy;
        }
 
        /// <summary>
        /// Finds all roots of a monic complex polynomial (highest degree first)
        /// using the Durand-Kerner (Weierstrass) iteration.
        /// </summary>
        private Complex[] FindRoots(Complex[] poly)
        {
            int n = poly.Length - 1;                 // degree
            var roots = new Complex[n];
 
            // Spread the initial guesses around a spiral to aid convergence.
            var seed = new Complex(0.4, 0.9);
            Complex cur = Complex.One;
            for (int i = 0; i < n; i++) { cur *= seed; roots[i] = cur; }
 
            for (int iter = 0; iter < 100; iter++)
            {
                double maxDelta = 0.0;
                for (int i = 0; i < n; i++)
                {
                    Complex num = EvalPoly(poly, roots[i]);
                    Complex den = Complex.One;
                    for (int j = 0; j < n; j++)
                        if (j != i) den *= (roots[i] - roots[j]);
 
                    if (den == Complex.Zero) continue;
                    Complex delta = num / den;
                    roots[i] -= delta;
                    double m = delta.Magnitude;
                    if (m > maxDelta) maxDelta = m;
                }
                if (maxDelta < 1e-12) break;          // converged
            }
            return roots;
        }
 
        /// <summary>Evaluates a polynomial (highest degree first) via Horner.</summary>
        private Complex EvalPoly(Complex[] poly, Complex z)
        {
            Complex result = Complex.Zero;
            for (int i = 0; i < poly.Length; i++)
                result = result * z + poly[i];
            return result;
        }
 
        /// <summary>
        /// Source-filter synthesis. A continuous excitation (impulse train for
        /// voiced frames, white noise for unvoiced) is run through the all-pole
        /// LPC filter. Output is de-emphasized and peak-normalized before writing.
        /// </summary>
        public void Synthesize(List<FrameInfo> frameData, int samplesPerFrame,
            int sampleRate, string fileName, bool useFixedPitch = true)
        {
            var rand = new Random(1234);
            int maxOrder = Order;
 
            // Filter history (continuous across frames -> no clicks at boundaries).
            var bp = new float[maxOrder];
            int offset = 0;
 
            int pulseCountdown = 0;   // samples until next glottal pulse
            var output = new List<float>();
            
            foreach (var fi in frameData)
            {
                float[] co = fi.Lpc;
                for (int smp = 0; smp < samplesPerFrame; smp++)
                {
                    // --- Excitation ---
                    float e;
                    if (fi.Voiced)
                    {
                        int T;
                        if (useFixedPitch)
                        {
                            T = FixedPitchHz > 0 ? sampleRate / FixedPitchHz : Math.Max(1, fi.PitchPeriod);
                        }
                        else
                        {
                            T = Math.Max(1, fi.PitchPeriod);
                            T -= PitchModulation;
                        }
                        if (pulseCountdown <= 0)
                        {
                            // Impulse scaled so average power matches the gain.
                            e = fi.Gain * (float)Math.Sqrt(T);
                            pulseCountdown = T;
                        }
                        else e = 0f;
                        pulseCountdown--;
                    }
                    else
                    {
                        e = fi.Gain * NextGaussian(rand);
                    }
 
                    // --- All-pole filter: y[n] = e - sum a_k * y[n-k] ---
                    float sum = e;
                    for (int j = 0; j < maxOrder; j++)
                    {
                        int indx = (offset + maxOrder - j) % maxOrder;
                        sum -= co[j] * bp[indx];
                    }
 
                    offset = (offset + 1) % maxOrder;
                    bp[offset] = sum;
                    output.Add(sum);
                }
            }
 
            // De-emphasis (inverse of the analysis pre-emphasis): y[n] = x[n] + a*y[n-1]
            float prev = 0f;
            for (int i = 0; i < output.Count; i++)
            {
                float v = output[i] + PreEmphasis * prev;
                output[i] = v;
                prev = v;
            }
 
            // Peak-normalize to avoid clipping.
            float peak = 0f;
            foreach (var v in output) peak = Math.Max(peak, Math.Abs(v));
            float scale = peak > 0f ? 0.95f / peak : 1f;
 
            var waveFormat = new WaveFormat(sampleRate, 16, 1); // mono, matches source rate
            using var writer = new WaveFileWriter(fileName, waveFormat);
            foreach (var v in output)
                writer.WriteSample(v * scale);
        }
 
        /// <summary>Standard normal sample via Box-Muller.</summary>
        private static float NextGaussian(Random rand)
        {
            double u1 = 1.0 - rand.NextDouble();
            double u2 = 1.0 - rand.NextDouble();
            return (float)(Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2));
        }
    }
}