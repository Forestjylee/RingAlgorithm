import { FilterType, Filter, Butterworth } from "./bfilter.js"

var samplingRate = 1000;
var filterOrder = 2;

var filt = new Butterworth(samplingRate, FilterType.Lowpass, filterOrder, [100]);
console.log(filt.Coefficients);

var filt = new Butterworth(samplingRate, FilterType.Highpass, filterOrder, [100]);
console.log(filt.Coefficients);

var filt = new Butterworth(samplingRate, FilterType.Bandpass, filterOrder, [100, 200]);
console.log(filt.Coefficients);

var filt = new Butterworth(samplingRate, FilterType.Notch, filterOrder, [48, 52]);
console.log(filt.Coefficients);