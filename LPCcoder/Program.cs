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

        static void HanningWindow()
        {
            float[] signalArray = new float[512];
            for (int i = 0; i < signalArray.Length; i++)
            {
                float hanningValue = (float)(0.5 * (1 - Math.Cos((2 * Math.PI * i) / 512)));

                Console.WriteLine($"{i}:\t {hanningValue.ToString()}");
                signalArray[i] *= hanningValue;
            }
        }

        static void Pinchettes()
        {
            float[] result =     [0.1f, 0.3f, 0.5f, 0.7f, 0.9f, 0.1f, 0.3f, 0.5f, 0.7f, 0.9f, 0.1f, 0.2f];
            float[] coffs = [0.01f, -0.02f, -0.03f, 0.04f, 0.01f, 0.06f, -0.007f, 0.01f, 0.02f, 0.01f, -0.02f, 0.03f];
            float[] buff =    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

            // e[10] is sound input and _s is buffer
            for (int x = 9; x >= 0; x--)
            {
                result[x] = result[x+1] - coffs[x] * buff[x];
            }
            for(int x = 8; x >= 0; x--)
            {
                buff[x+1] = buff[x] + coffs[x] * result[x];
            }
            buff[0] = result[0];
            for(int y = 9; y >= 0; y--)
            {
                Console.WriteLine($"{y} : e:{result[y]}\t\t s_{buff[y]}");
            }

            //*excitation++ = e[10];
            //*output++ = e[0];
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
        public (List<float[]>,List<float>) CreateFrames(float[]input, int frameSize)
        {
            int sampleLength = input.Length;
            int totalFrames = sampleLength / frameSize;
            List<float[]> result = [];  
            List<float> volume = [];
            int c = 0;
            float average = 0.0f;
            for (int y = 0; y < totalFrames; y++)
            {
                float[] frame = new float[frameSize];
                for (int x = 0; x < frameSize; x++)
                {
                    average += frame[x] = input[x+c];
                    frame[x] = frame[x] * CalculateHanning(x, frameSize);
                }
                c += frameSize;
                average /= frameSize; // average vol
                result.Add(frame);
                volume.Add(average);
            }
            return (result,volume);
        }

        /// <summary>
        /// räknar ut Hanningvärdet baserad på fönstrets längd och position
        /// </summary>
        /// <param name="index">position</param>
        /// <param name="length">fönstrets längd</param>
        /// <returns></returns>
        public float CalculateHanning(int index, int length)
        {
            float hanningValue = (float)(0.5 * (1 - Math.Cos((2 * Math.PI * index) / length)));
            return hanningValue;
        }

        /// <summary>
        /// Cuts all zero:es at the beginning of the soundfile
        /// </summary>
        /// <param name="input"></param>
        /// <returns></returns>
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
            int VoicePitch = 1;
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
                        //// Generate buzz, at the specified voice pitch
                        count++;
                        float w = (float)(count %= rate / VoicePitch) / (rate / VoicePitch);

                        float pt = (float)Math.Pow(2.0, w);
                        float f = (float)(pt - 1 / (1 + w)); // -0.5  + rand.NextDouble()) * 0.5  + pt - 1 / (1 + w)

                        // Apply the filter (LPC coefficients) co[j]
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
            List<float[]> coffs = [];
            
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
            Synthesize(coffs, 24, order);
            
            Console.WriteLine("Your soundfile is ready");
            
            
        }
    }
}
