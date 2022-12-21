/**
This library is based on mkfilter.c by Anthony J. Fisher, University of York, September 1992.
see:
http://www-users.cs.york.ac.uk/~fisher/mkfilter
https://github.com/university-of-york/cs-www-users-fisher
*/
var __classPrivateFieldGet = (this && this.__classPrivateFieldGet) || function (receiver, state, kind, f) {
    if (kind === "a" && !f) throw new TypeError("Private accessor was defined without a getter");
    if (typeof state === "function" ? receiver !== state || !f : !state.has(receiver)) throw new TypeError("Cannot read private member from an object whose class did not declare it");
    return kind === "m" ? f : kind === "a" ? f.call(receiver) : f ? f.value : state.get(receiver);
};
var __classPrivateFieldSet = (this && this.__classPrivateFieldSet) || function (receiver, state, value, kind, f) {
    if (kind === "m") throw new TypeError("Private method is not writable");
    if (kind === "a" && !f) throw new TypeError("Private accessor was defined without a setter");
    if (typeof state === "function" ? receiver !== state || !f : !state.has(receiver)) throw new TypeError("Cannot write private member to an object whose class did not declare it");
    return (kind === "a" ? f.call(receiver, value) : f ? f.value = value : state.set(receiver, value)), value;
};
var _Butterworth_instances, _Butterworth_Cmone, _Butterworth_Czero, _Butterworth_Cone, _Butterworth_Ctwo, _Butterworth_Chalf, _Butterworth_expand, _Butterworth_multin, _Butterworth_evaluate, _Butterworth_eval, _Butterworth_xsqrt, _Butterworth_csqrt, _Butterworth_cexp, _Butterworth_cadd, _Butterworth_csub, _Butterworth_cmul, _Butterworth_cdiv, _Butterworth_cneg, _a, _Filter_filt, _RealtimeFilter_coeff, _RealtimeFilter_channels, _RealtimeFilter_x, _RealtimeFilter_y;
/**
 * Enumeration containing available filter types.
 */
export var FilterType;
(function (FilterType) {
    FilterType[FilterType["Lowpass"] = 0] = "Lowpass";
    FilterType[FilterType["Highpass"] = 1] = "Highpass";
    FilterType[FilterType["Bandpass"] = 2] = "Bandpass";
    FilterType[FilterType["Notch"] = 3] = "Notch";
})(FilterType || (FilterType = {}));
/**
 * Filter coefficients
 */
export class FilterCoefficients {
    /**
     * Initializes a new instance of @see FilterCoefficient.
     * @param a a coefficients.
     * @param b b coefficients.
     */
    constructor(a, b) {
        this.a = a;
        this.b = b;
    }
}
/**
 * Class to calculate Butterworth IIR filter coefficients.
 */
export class Butterworth {
    //#endregion
    /**
     * Initializes a new instance of @see Butterworth.
     * @param samplingRate The sampling rate of the signal.
     * @param type Parameter to select the filter behavior @see FilterType.
     * @param order The filter order.
     * @param cutoffFrequencies An array with one filter cutoff frequency for @see FilterType.Highpass and @see FilterType.Lowpass. An array with two filter cutoff frequencies for @type FilterType.Bandpass and @type FilterType.Notch.
     */
    constructor(samplingRate, type, order, cutoffFrequencies) {
        _Butterworth_instances.add(this);
        //#region  private members...
        _Butterworth_Cmone.set(this, new Complex(-1.0, 0.0));
        _Butterworth_Czero.set(this, new Complex(0.0, 0.0));
        _Butterworth_Cone.set(this, new Complex(1.0, 0.0));
        _Butterworth_Ctwo.set(this, new Complex(2.0, 0.0));
        _Butterworth_Chalf.set(this, new Complex(0.5, 0.0));
        var cutoffLow = -1;
        var cutoffHigh = -1;
        if ((type == FilterType.Lowpass || type == FilterType.Highpass) && cutoffFrequencies.length != 1)
            throw new Error("Only one cutoff frequency allowed for low- and highpass filters.");
        if ((type == FilterType.Bandpass || type == FilterType.Notch) && cutoffFrequencies.length != 2)
            throw new Error("Two cutoff frequencies required for bandpass and notch filters.");
        var cutoffLow = cutoffFrequencies[0];
        if (type == FilterType.Bandpass || type == FilterType.Notch)
            var cutoffHigh = cutoffFrequencies[1];
        if ((type == FilterType.Bandpass || type == FilterType.Notch) && cutoffLow >= cutoffHigh)
            throw new Error("Lower cutoff frequency mut be lower than higher cutoff frequency");
        if ((type == FilterType.Lowpass || type == FilterType.Highpass) && cutoffHigh != -1)
            throw new Error("Only lower cutoff frequency must be set for low- or highpass filter types.");
        if ((type == FilterType.Lowpass || type == FilterType.Highpass) && cutoffLow > samplingRate / 2)
            throw new Error("Cutoff must be lower than half of the sampling rate.");
        if ((type == FilterType.Bandpass || type == FilterType.Notch) && (cutoffLow > samplingRate / 2 || cutoffHigh > samplingRate / 2))
            throw new Error("Cutoff must be lower than half of the sampling rate.");
        //set defaults
        var numpoles = 0;
        var size = 0;
        if (type == FilterType.Bandpass || type == FilterType.Notch)
            size = 2 * order;
        else
            size = order;
        var spoles = new Array(size);
        var zpoles = new Array(size);
        var zzeros = new Array(size);
        //compute s
        for (var i = 0; i < 2 * order; i++) {
            var s = new Complex(0.0, (order & 1) == 1 ? (i * Math.PI) / order : ((i + 0.5) * Math.PI) / order);
            var z = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cexp).call(this, s);
            if (z.Real < 0.0)
                spoles[numpoles++] = z;
        }
        //normalize
        var rawAlpha1 = 0.0;
        var rawAlpha2 = 0.0;
        rawAlpha1 = cutoffLow / samplingRate;
        if (type == FilterType.Bandpass || type == FilterType.Notch)
            rawAlpha2 = cutoffHigh / samplingRate;
        else
            rawAlpha2 = rawAlpha1;
        var warpedAlpha1 = Math.tan(Math.PI * rawAlpha1) / Math.PI;
        var warpedAlpha2 = Math.tan(Math.PI * rawAlpha2) / Math.PI;
        var w1 = new Complex(2.0 * Math.PI * warpedAlpha1, 0.0);
        var w2 = new Complex(2.0 * Math.PI * warpedAlpha2, 0.0);
        var w0 = __classPrivateFieldGet(this, _Butterworth_Czero, "f");
        var bw = __classPrivateFieldGet(this, _Butterworth_Czero, "f");
        switch (type) {
            case FilterType.Lowpass:
                for (var i = 0; i < numpoles; i++)
                    spoles[i] = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cmul).call(this, spoles[i], w1);
                break;
            case FilterType.Highpass:
                for (var i = 0; i < numpoles; i++)
                    spoles[i] = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cdiv).call(this, w1, spoles[i]);
                break;
            case FilterType.Bandpass:
                w0 = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_csqrt).call(this, __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cmul).call(this, w1, w2));
                bw = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_csub).call(this, w2, w1);
                for (var i = 0; i < numpoles; i++) {
                    var hba = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cmul).call(this, __classPrivateFieldGet(this, _Butterworth_Chalf, "f"), __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cmul).call(this, spoles[i], bw));
                    var temp = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cdiv).call(this, w0, hba);
                    temp = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_csqrt).call(this, __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_csub).call(this, __classPrivateFieldGet(this, _Butterworth_Cone, "f"), __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cmul).call(this, temp, temp)));
                    spoles[i] = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cmul).call(this, hba, __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cadd).call(this, __classPrivateFieldGet(this, _Butterworth_Cone, "f"), temp));
                    spoles[numpoles + i] = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cmul).call(this, hba, __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_csub).call(this, __classPrivateFieldGet(this, _Butterworth_Cone, "f"), temp));
                }
                numpoles *= 2;
                break;
            case FilterType.Notch:
                w0 = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_csqrt).call(this, __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cmul).call(this, w1, w2));
                bw = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_csub).call(this, w2, w1);
                for (var i = 0; i < numpoles; i++) {
                    var hba = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cmul).call(this, __classPrivateFieldGet(this, _Butterworth_Chalf, "f"), __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cdiv).call(this, bw, spoles[i]));
                    var temp = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cdiv).call(this, w0, hba);
                    temp = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_csqrt).call(this, __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_csub).call(this, __classPrivateFieldGet(this, _Butterworth_Cone, "f"), __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cmul).call(this, temp, temp)));
                    spoles[i] = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cmul).call(this, hba, __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cadd).call(this, __classPrivateFieldGet(this, _Butterworth_Cone, "f"), temp));
                    spoles[numpoles + i] = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cmul).call(this, hba, __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_csub).call(this, __classPrivateFieldGet(this, _Butterworth_Cone, "f"), temp));
                }
                numpoles *= 2;
                break;
            default:
                throw new Error("Unknown filter type");
        }
        //compute z
        for (var i = 0; i < numpoles; i++) {
            var top = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cadd).call(this, __classPrivateFieldGet(this, _Butterworth_Ctwo, "f"), spoles[i]);
            var bot = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_csub).call(this, __classPrivateFieldGet(this, _Butterworth_Ctwo, "f"), spoles[i]);
            zpoles[i] = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cdiv).call(this, top, bot);
            switch (type) {
                case FilterType.Lowpass:
                    zzeros[i] = __classPrivateFieldGet(this, _Butterworth_Cmone, "f");
                    break;
                case FilterType.Highpass:
                    zzeros[i] = __classPrivateFieldGet(this, _Butterworth_Cone, "f");
                    break;
                case FilterType.Bandpass:
                    if (i % 2 == 0)
                        zzeros[i] = __classPrivateFieldGet(this, _Butterworth_Cone, "f");
                    else
                        zzeros[i] = __classPrivateFieldGet(this, _Butterworth_Cmone, "f");
                    break;
                case FilterType.Notch:
                    if (i % 2 == 0)
                        zzeros[i] = new Complex(Math.cos(2 * Math.PI * ((cutoffHigh + cutoffLow) / 2) / samplingRate), Math.sin(2 * Math.PI * ((cutoffHigh + cutoffLow) / 2) / samplingRate));
                    else
                        zzeros[i] = new Complex(Math.cos(2 * Math.PI * ((cutoffHigh + cutoffLow) / 2) / samplingRate), -Math.sin(2 * Math.PI * ((cutoffHigh + cutoffLow) / 2) / samplingRate));
                    break;
                default:
                    throw new Error("Unknown filter type");
            }
        }
        //expand poly
        var topCoeffs = new Array(numpoles + 1);
        var botCoeffs = new Array(numpoles + 1);
        var a = new Array(numpoles + 1);
        var b = new Array(numpoles + 1);
        __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_expand).call(this, zzeros, topCoeffs, numpoles);
        __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_expand).call(this, zpoles, botCoeffs, numpoles);
        var FCgain = __classPrivateFieldGet(this, _Butterworth_Cone, "f");
        if (type != FilterType.Notch) {
            var st = new Complex(0.0, Math.PI * (rawAlpha1 + rawAlpha2));
            var zfc = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cexp).call(this, st);
            FCgain = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_evaluate).call(this, topCoeffs, botCoeffs, numpoles, zfc);
        }
        else
            FCgain = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_evaluate).call(this, topCoeffs, botCoeffs, numpoles, __classPrivateFieldGet(this, _Butterworth_Cone, "f"));
        for (var i = 0; i <= numpoles; i++) {
            if (type == FilterType.Highpass || type == FilterType.Lowpass)
                b[i] = topCoeffs[i].Real / botCoeffs[numpoles].Real / FCgain.Magnitude / Math.sqrt(2);
            else
                b[i] = topCoeffs[i].Real / botCoeffs[numpoles].Real / FCgain.Magnitude;
            a[i] = botCoeffs[i].Real / botCoeffs[numpoles].Real;
        }
        this.Coefficients = new FilterCoefficients(a.reverse(), b.reverse());
    }
}
_Butterworth_Cmone = new WeakMap(), _Butterworth_Czero = new WeakMap(), _Butterworth_Cone = new WeakMap(), _Butterworth_Ctwo = new WeakMap(), _Butterworth_Chalf = new WeakMap(), _Butterworth_instances = new WeakSet(), _Butterworth_expand = function _Butterworth_expand(pz, coeffs, numpoles) {
    coeffs[0] = __classPrivateFieldGet(this, _Butterworth_Cone, "f");
    for (var i = 0; i < numpoles; i++)
        coeffs[i + 1] = __classPrivateFieldGet(this, _Butterworth_Czero, "f");
    for (var i = 0; i < numpoles; i++)
        __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_multin).call(this, pz[i], coeffs, numpoles);
    for (var i = 0; i < numpoles + 1; i++) {
        if (Math.abs(coeffs[i].Imaginary) > 1e-10)
            throw new Error("Filter calculation failed. Coefficients of z^k not real.");
    }
}, _Butterworth_multin = function _Butterworth_multin(w, coeffs, numpoles) {
    var nw = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cneg).call(this, w);
    for (var i = numpoles; i >= 1; i--)
        coeffs[i] = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cadd).call(this, __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cmul).call(this, nw, coeffs[i]), coeffs[i - 1]);
    coeffs[0] = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cmul).call(this, nw, coeffs[0]);
}, _Butterworth_evaluate = function _Butterworth_evaluate(topco, botco, np, z) {
    return __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cdiv).call(this, __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_eval).call(this, topco, np, z), __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_eval).call(this, botco, np, z));
}, _Butterworth_eval = function _Butterworth_eval(coeffs, np, z) {
    var sum = __classPrivateFieldGet(this, _Butterworth_Czero, "f");
    for (var i = np; i >= 0; i--)
        sum = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cadd).call(this, __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_cmul).call(this, sum, z), coeffs[i]);
    return sum;
}, _Butterworth_xsqrt = function _Butterworth_xsqrt(x) {
    return (x >= 0.0) ? Math.sqrt(x) : 0.0;
}, _Butterworth_csqrt = function _Butterworth_csqrt(x) {
    var r = Math.sqrt(Math.pow(x.Real, 2) + Math.pow(x.Imaginary, 2));
    var real = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_xsqrt).call(this, 0.5 * (r + x.Real));
    var imag = __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_xsqrt).call(this, 0.5 * (r - x.Real));
    if (x.Imaginary < 0.0)
        imag = -imag;
    return new Complex(real, imag);
}, _Butterworth_cexp = function _Butterworth_cexp(x) {
    var r = Math.exp(x.Real);
    return new Complex(r * Math.cos(x.Imaginary), r * Math.sin(x.Imaginary));
}, _Butterworth_cadd = function _Butterworth_cadd(x, y) {
    return new Complex(x.Real + y.Real, x.Imaginary + y.Imaginary);
}, _Butterworth_csub = function _Butterworth_csub(x, y) {
    return new Complex(x.Real - y.Real, x.Imaginary - y.Imaginary);
}, _Butterworth_cmul = function _Butterworth_cmul(x, y) {
    return new Complex((x.Real * y.Real - x.Imaginary * y.Imaginary), (x.Imaginary * y.Real + x.Real * y.Imaginary));
}, _Butterworth_cdiv = function _Butterworth_cdiv(x, y) {
    var mag = y.Real * y.Real + y.Imaginary * y.Imaginary;
    return new Complex((x.Real * y.Real + x.Imaginary * y.Imaginary) / mag, (x.Imaginary * y.Real - x.Real * y.Imaginary) / mag);
}, _Butterworth_cneg = function _Butterworth_cneg(x) {
    return __classPrivateFieldGet(this, _Butterworth_instances, "m", _Butterworth_csub).call(this, __classPrivateFieldGet(this, _Butterworth_Czero, "f"), x);
};
export class Filter {
    /**
     * Applies a filter to the input data.
     * @param data 2D Array structured as [samples, channels].
     * @param filt A filter object @see Butterworth.
     * @returns Filtered data.
     */
    static filter(data, filt) {
        return __classPrivateFieldGet(this, _a, "m", _Filter_filt).call(this, data, filt, false, false);
    }
    /**
     * Applies a filter to the input data twice; forwards and backwards.
     * @param data 2D Array structured as [samples, channels].
     * @param filt A filter object @see Butterworth.
     * @returns Filtered data.
     */
    static filtfilt(data, filt) {
        return __classPrivateFieldGet(this, _a, "m", _Filter_filt).call(this, data, filt, true, false);
    }
}
_a = Filter, _Filter_filt = function _Filter_filt(data, filt, filtfilt, offsetCorrection) {
    var rows = data.length;
    var columns = data[0].length;
    var coeff = filt.Coefficients;
    //allocate buffers
    var dataOut = new Array(rows);
    for (var i = 0; i < rows; i++) {
        dataOut[i] = new Array(columns);
    }
    if (coeff.a.length != coeff.b.length)
        throw new Error("Invalid filter coefficients.");
    var numberOfCoefficients = coeff.a.length;
    var x = new Array(numberOfCoefficients);
    var y = new Array(numberOfCoefficients);
    for (var i = 0; i < numberOfCoefficients; i++) {
        x[i] = new Array(columns);
    }
    for (var i = 0; i < numberOfCoefficients; i++) {
        y[i] = new Array(columns);
    }
    var xyr = x.length;
    var xyc = x[0].length;
    for (var i = 0; i < xyr; i++) {
        for (var j = 0; j < xyc; j++) {
            if (offsetCorrection) {
                x[i][j] = data[0][0];
                y[i][j] = data[0][0];
            }
            else {
                x[i][j] = 0;
                y[i][j] = 0;
            }
        }
    }
    //filter
    for (var c = 0; c < columns; c++) {
        for (var r = 0; r < rows; r++) {
            //shift buffer
            for (var i = 0; i < numberOfCoefficients - 1; i++) {
                x[i][c] = x[i + 1][c];
                y[i][c] = y[i + 1][c];
            }
            //transfer function
            x[numberOfCoefficients - 1][c] = data[r][c];
            y[numberOfCoefficients - 1][c] = coeff.b[0] * x[numberOfCoefficients - 1][c];
            for (var i = 1; i < numberOfCoefficients; i++) {
                y[numberOfCoefficients - 1][c] = y[numberOfCoefficients - 1][c] + coeff.b[i] * x[numberOfCoefficients - 1 - i][c] - coeff.a[i] * y[numberOfCoefficients - 1 - i][c];
            }
            dataOut[r][c] = y[numberOfCoefficients - 1][c];
        }
    }
    //filter reverse
    if (filtfilt) {
        for (var c = columns - 1; c >= 0; c--) {
            for (var r = rows - 1; r >= 0; r--) {
                //shift buffer
                for (var i = 0; i < numberOfCoefficients - 1; i++) {
                    x[i][c] = x[i + 1][c];
                    y[i][c] = y[i + 1][c];
                }
                //transfer function
                x[numberOfCoefficients - 1][c] = dataOut[r][c];
                y[numberOfCoefficients - 1][c] = coeff.b[0] * x[numberOfCoefficients - 1][c];
                for (var i = 1; i < numberOfCoefficients; i++) {
                    y[numberOfCoefficients - 1][c] = y[numberOfCoefficients - 1][c] + coeff.b[i] * x[numberOfCoefficients - 1 - i][c] - coeff.a[i] * y[numberOfCoefficients - 1 - i][c];
                }
                dataOut[r][c] = y[numberOfCoefficients - 1][c];
            }
        }
    }
    return dataOut;
};
export class RealtimeFilter {
    constructor(filt, channels) {
        _RealtimeFilter_coeff.set(this, void 0);
        _RealtimeFilter_channels.set(this, void 0);
        _RealtimeFilter_x.set(this, void 0);
        _RealtimeFilter_y.set(this, void 0);
        __classPrivateFieldSet(this, _RealtimeFilter_coeff, filt.Coefficients, "f");
        __classPrivateFieldSet(this, _RealtimeFilter_channels, channels, "f");
        if (__classPrivateFieldGet(this, _RealtimeFilter_coeff, "f").a.length != __classPrivateFieldGet(this, _RealtimeFilter_coeff, "f").b.length)
            throw new Error("Invalid filter coefficients.");
        var numberOfCoefficients = __classPrivateFieldGet(this, _RealtimeFilter_coeff, "f").a.length;
        __classPrivateFieldSet(this, _RealtimeFilter_x, new Array(numberOfCoefficients), "f");
        __classPrivateFieldSet(this, _RealtimeFilter_y, new Array(numberOfCoefficients), "f");
        for (var i = 0; i < numberOfCoefficients; i++) {
            __classPrivateFieldGet(this, _RealtimeFilter_x, "f")[i] = new Array(__classPrivateFieldGet(this, _RealtimeFilter_channels, "f"));
        }
        for (var i = 0; i < numberOfCoefficients; i++) {
            __classPrivateFieldGet(this, _RealtimeFilter_y, "f")[i] = new Array(__classPrivateFieldGet(this, _RealtimeFilter_channels, "f"));
        }
        var xyr = __classPrivateFieldGet(this, _RealtimeFilter_x, "f").length;
        var xyc = __classPrivateFieldGet(this, _RealtimeFilter_x, "f")[0].length;
        for (var i = 0; i < xyr; i++) {
            for (var j = 0; j < xyc; j++) {
                __classPrivateFieldGet(this, _RealtimeFilter_x, "f")[i][j] = 0;
                __classPrivateFieldGet(this, _RealtimeFilter_y, "f")[i][j] = 0;
            }
        }
    }
    step(data) {
        var rows = data.length;
        var columns = data[0].length;
        var numberOfCoefficients = __classPrivateFieldGet(this, _RealtimeFilter_coeff, "f").a.length;
        if (columns != __classPrivateFieldGet(this, _RealtimeFilter_channels, "f"))
            throw new Error('Invalid data dimensions.');
        var dataOut = new Array(rows);
        for (var i = 0; i < rows; i++) {
            dataOut[i] = new Array(columns);
        }
        for (var c = 0; c < columns; c++) {
            for (var r = 0; r < rows; r++) {
                //shift buffer
                for (var i = 0; i < numberOfCoefficients - 1; i++) {
                    __classPrivateFieldGet(this, _RealtimeFilter_x, "f")[i][c] = __classPrivateFieldGet(this, _RealtimeFilter_x, "f")[i + 1][c];
                    __classPrivateFieldGet(this, _RealtimeFilter_y, "f")[i][c] = __classPrivateFieldGet(this, _RealtimeFilter_y, "f")[i + 1][c];
                }
                //transfer function
                __classPrivateFieldGet(this, _RealtimeFilter_x, "f")[numberOfCoefficients - 1][c] = data[r][c];
                __classPrivateFieldGet(this, _RealtimeFilter_y, "f")[numberOfCoefficients - 1][c] = __classPrivateFieldGet(this, _RealtimeFilter_coeff, "f").b[0] * __classPrivateFieldGet(this, _RealtimeFilter_x, "f")[numberOfCoefficients - 1][c];
                for (var i = 1; i < numberOfCoefficients; i++) {
                    __classPrivateFieldGet(this, _RealtimeFilter_y, "f")[numberOfCoefficients - 1][c] = __classPrivateFieldGet(this, _RealtimeFilter_y, "f")[numberOfCoefficients - 1][c] + __classPrivateFieldGet(this, _RealtimeFilter_coeff, "f").b[i] * __classPrivateFieldGet(this, _RealtimeFilter_x, "f")[numberOfCoefficients - 1 - i][c] - __classPrivateFieldGet(this, _RealtimeFilter_coeff, "f").a[i] * __classPrivateFieldGet(this, _RealtimeFilter_y, "f")[numberOfCoefficients - 1 - i][c];
                }
                dataOut[r][c] = __classPrivateFieldGet(this, _RealtimeFilter_y, "f")[numberOfCoefficients - 1][c];
            }
        }
        return dataOut;
    }
}
_RealtimeFilter_coeff = new WeakMap(), _RealtimeFilter_channels = new WeakMap(), _RealtimeFilter_x = new WeakMap(), _RealtimeFilter_y = new WeakMap();
class Complex {
    constructor(real, imaginary) {
        this.Real = real;
        this.Imaginary = imaginary;
        this.Magnitude = Math.sqrt(Math.pow(this.Real, 2) + Math.pow(this.Imaginary, 2));
    }
}