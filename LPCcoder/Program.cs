namespace LPCcoder
{
    using Microsoft.VisualBasic;
    using NAudio.Wave;
    using System.Linq.Expressions;
    using static System.Runtime.InteropServices.JavaScript.JSType;

    internal class Program
    {
        
        static void Main(string[] args)
        {
            LPC lpc = new();
            lpc.StartHere("testmening3.wav");
        }
    }
    
    public class LPC
    {
        const string FileName = "Result.wav";
        public float[] ReadWav(string filename)
        {
            AudioFileReader reader = new AudioFileReader(filename);
            ISampleProvider isp = reader.ToSampleProvider();
            float[] buffer = new float[reader.Length / 2];
            isp.Read(buffer, 0, buffer.Length);
            return buffer;
        }

        /// <summary>
        /// Old method
        /// </summary>
        /// <param name="signal"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        /// 
        [Obsolete]
        public static float[] CalculateLPC(float[] signal, int order)
        {
            int n = signal.Length;
            float[] autocorrelation = new float[order + 1];
            float[] lpcCoefficients = new float[order + 1];
            float[] error = new float[order + 1];
            float[] reflection = new float[order + 1];

            // Calculate autocorrelation
            for (int i = 0; i <= order; i++)
            {
                autocorrelation[i] = 0.0f;
                for (int j = 0; j < n - i; j++)
                {
                    autocorrelation[i] += signal[j] * signal[j + i];
                }
            }

            // Initialize error
            error[0] = autocorrelation[0];

            // Levinson-Durbin recursion
            for (int i = 1; i <= order; i++)
            {
                float sum = 0.0f;
                for (int j = 1; j < i; j++)
                {
                    sum += lpcCoefficients[j] * autocorrelation[i - j];
                }
                reflection[i] = (autocorrelation[i] - sum) / error[i - 1];
                lpcCoefficients[i] = reflection[i];

                for (int j = 1; j < i; j++)
                {
                    lpcCoefficients[j] -= reflection[i] * lpcCoefficients[i - j];
                }

                error[i] = (1.0f - reflection[i] * reflection[i]) * error[i - 1];
            }

            return lpcCoefficients;
        }

        public static void PerformLPCAnalysis(float[] frame, int order, float[] lpcCoeffs, float[] autocorr, float[] reflectionCoeffs)
        {
            var frameSize = frame.Length;
            // Step 1: Calculate the autocorrelation coefficients
            for (int k = 0; k <= order; k++)
            {
                float sum = 0.0f;
                for (int n = 0; n < frameSize - k; n++)
                {
                    sum += frame[n] * frame[n + k];
                }
                autocorr[k] = sum;
            }

            // Step 2: Solve the Toeplitz system of equations using Levinson-Durbin recursion
            lpcCoeffs[0] = 1.0f;
            float error = autocorr[0];

            if (error == 0.0f)
            {
                // If the autocorrelation is zero, return zero coefficients
                for (int i = 0; i <= order; i++)
                {
                    lpcCoeffs[i] = 0.0f;
                }
                return;
            }

            reflectionCoeffs[0] = -autocorr[1] / autocorr[0];
            lpcCoeffs[1] = reflectionCoeffs[0];
            error *= (1.0f - reflectionCoeffs[0] * reflectionCoeffs[0]);

            for (int m = 1; m < order; m++)
            {
                float sum = 0.0f;
                for (int j = 0; j <= m; j++)
                {
                    sum += lpcCoeffs[j] * autocorr[m + 1 - j];
                }

                reflectionCoeffs[m] = -sum / error;

                for (int j = 1; j <= (m + 1) / 2; j++)
                {
                    float tmp = lpcCoeffs[j] + reflectionCoeffs[m] * lpcCoeffs[m + 1 - j];
                    lpcCoeffs[m + 1 - j] += reflectionCoeffs[m] * lpcCoeffs[j];
                    lpcCoeffs[j] = tmp;
                }
                lpcCoeffs[m + 1] = reflectionCoeffs[m];
                error *= (1.0f - reflectionCoeffs[m] * reflectionCoeffs[m]);
            }

            // At this point, lpcCoeffs contains the LPC coefficients and
            // reflectionCoeffs contains the reflection coefficients.
        }
        public (List<float[]>,List<float>) CreateFrames(float[]input, int size)
        {
            int n = input.Length;
            int resultAmounts = n / size;
            List<float[]> result = new List<float[]>();  
            List<float> volume = new List<float>();
            int c = 0;
            float average = 0.0f;
            for (int y = 0; y < resultAmounts; y++)
            {
                float[] frame = new float[size];
                for (int x = 0; x < size; x++)
                {
                    average += frame[x] = input[x+c];
                }
                c += size;
                average /= size; // average vol
                result.Add(frame);
                volume.Add(average);
            }
            return (result,volume);
        }
        public float[] CleanUp(float[] input)
        {
            int index = 0;
            while (index < input.Length && input[index] == 0.0f) index++;
            return input.Skip(index).ToArray();
        }

        public void Synthesize(List<float[]> coffs, int lengthMultiplicator, int Order)
        {
            var rand = new Random();
            int sampleRate = 44100;
            int VoicePitch = 30;
            int count = 0;
            float[] k = new float[Order];
            float[] bp = new float[Order];
            int offset = 0;
            int MaxOrder = Order;
            //float phase = 0;
            //float fTablesize = 512;
            float sampling_period = 0.00002267573696f;
            WaveFormat waveFormat = new WaveFormat();
            using (WaveFileWriter writer = new WaveFileWriter(FileName, waveFormat))
            {
                foreach (var co in coffs)
                {
                    int rate = (int)(1.0f / sampling_period); 
                    int samples_per_frame = (int)(0.005 / sampling_period);
                    for (int smp = 0; smp < samples_per_frame * lengthMultiplicator; smp++)
                    {
                        // https://github.com/bisqwit/speech_synth_series/blob/master/ep4-speechsyn/pcmaudio-lpc-wav.cc

                        //// Generate buzz, at the specified voice pitch
                        count++;
                        float w = (float)(count %= rate / VoicePitch) / (rate / VoicePitch);

                        float pt = (float)Math.Pow(2.0, w);
                        float f = (float)(pt - 1 / (1 + w)); // -0.5  + rand.NextDouble()) * 0.5  + pt - 1 / (1 + w)

                        // Apply the filter (LPC coefficients)
                        float sum = f;
                        for (int j = 0; j < Order; j++)
                        {
                            int indx = (offset + MaxOrder - j) % MaxOrder;
                            sum -= (co[j]) * bp[indx];
                        }

                        // Save it into a rolling buffer
                        offset++;
                        int index = offset %= MaxOrder;
                        float r = bp[index] = sum;
                        writer.WriteSample(r);
                    }
                    
                }
            }
        }

        public void StartHere(string filePath)
        {
            int FrameLength = 2048;
            float[] rawSignal = ReadWav(filePath);
            float[] signal = CleanUp(rawSignal);
            int order = 20; // LPC order
            
            //float[] lpcCoefficients = CalculateLPC(signal, order);
            var (frames,volume) = CreateFrames(signal, FrameLength);
            List<float[]> coffs = new();
            foreach (var fr in frames)
            {
                float[] lpcCoeffs = new float[order+1];
                float[] autocorr = new float[order + 1];
                float[] reflectionCoeffs = new float[order + 1];
                PerformLPCAnalysis(fr, order, lpcCoeffs, autocorr, reflectionCoeffs);
                var lpcCoeffs2 = lpcCoeffs.Skip(1).
                    Take(order).
                    ToArray();
                coffs.Add(lpcCoeffs2);
            }
            Synthesize(coffs, 4, order);
            
            Console.WriteLine("Your synthesized speech is ready");
            
            
        }
    }
}
