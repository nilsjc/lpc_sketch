# LPC experiments with C# on .NET
![Alt text](lpcinuse.png?raw=true "Schematic")<br>
This C# console program will generate synthetic speech, based on a soundfile with recorded speech.<br>
When you start, it will ask for a file and some user input like formant or pitch shifting, or if you want your speech sound like a robot/vocoder. It will then generate a result.wav file<br>
If you are not satisifed, you can always try again. You will load the .wav-file only once, at the start of the program, and then you can retry the generator with different parameters. You will find the sound in the generated "result.wav"file.<br>
For now this sound quite Lofi, but i am anyway more interested in using this in a musical application.<br>
This is an work in progress and some features is not implemented, like time stretch, external source modulation or control.<br>
Both the pitch change and the formant shift is a bit shaky as well. Sometimes you will get only silence. I  recommend starting with all settings at default; that means no robot voice, voice/unvoiced mode, pitch = 0 and formant shift = 1.0. Then you can try to change the parameters to more extreme settings.<br>
I will improve this program in the near future. A GUI would also be nice to have, with some sound playing capabilities. I am using Naudio library that is quite old today (and not maintained to my knowing), however it works, but only on Windows i think. I will soon replace Naudio with Soundflow instead.<br>
In the not so near future this will be converted to a MCU, maybe an 32bit ARM platform, for use in a modular synthesizer. The plan is then to process speech in realtime, not rendered like now.<br>
Until then, this project is my test platform.<br>
Have fun..