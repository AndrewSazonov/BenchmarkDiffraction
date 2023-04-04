import sys, os
import timeit
import numpy as np
import scipy as sc
import lmfit
import iminuit

from PySide6.QtCore import QObject, Signal, Slot, Property
from PySide6.QtWidgets import QApplication
from PySide6.QtQml import QQmlApplicationEngine
from PySide6.QtGui import QSurfaceFormat
from PySide6.QtQuick3D import QQuick3D
from PySide6 import QtCharts

import CFML_api.powder_mod
import cryspy
from cryspy.procedure_rhochi.rhochi_by_dictionary import \
    rhochi_calc_chi_sq_by_dictionary


DTYPE = np.float64
RCIF_FILE_NAME = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'examples', 'PbSO4.rcif')


class Proxy(QObject):
    axisRangesChanged = Signal()
    timeChanged = Signal()
    measSerieChanged = Signal()
    calcSerieChanged = Signal()
    residSerieChanged = Signal()
    dictDataChanged = Signal(bool)
    wavelengthChanged = Signal(bool)
    resolutionXChanged = Signal(bool)

    meas2thetaOffsetChanged = Signal(bool)
    scaleChanged = Signal(bool)

    stepChanged = Signal(bool)
    arraySizeChanged = Signal()
    radiationProbeChanged = Signal(bool)
    calculatorChanged = Signal(bool)

    def __init__(self, parent=None):
        super(Proxy, self).__init__(parent)

        # CrysPy objects
        self._cryspyObj = cryspy.load_file(RCIF_FILE_NAME)
        self._cryspyDict = self._cryspyObj.get_dictionary()
        self._cryspyInOutDict = {}

        # EasyDiffraction objects
        self._edDict = self.createEdDict(self._cryspyObj, self._cryspyDict)

        # CrysFML objects
        self._crysfmlDict = self._edDict.copy()

        # Misc
        self._first_experiment_name = list(self._edDict['experiments'][0].keys())[0]

        self._wavelength = self._edDict['experiments'][0][self._first_experiment_name]['_diffrn_radiation_wavelength']
        self._resolutionX = self._edDict['experiments'][0][self._first_experiment_name]['_pd_instr_resolution_x']
        self._step = self._edDict['experiments'][0][self._first_experiment_name]['_pd_meas_2theta_range_inc']
        self._radiationProbe = self._edDict['experiments'][0][self._first_experiment_name]['_diffrn_radiation_probe']

        self._meas2thetaOffset = self._edDict['experiments'][0][self._first_experiment_name]['_pd_meas_2theta_offset']
        self._scale = self._edDict['experiments'][0][self._first_experiment_name]['_phase'][0]['_scale']

        self._calculator = 'cryspy'

        self._arraySize = 0
        self._axisRanges = {'xMin': 0, 'xMax': 1, 'yMin': 0, 'yMax': 1}
        self._measSerie = QtCharts.QXYSeries
        self._calcSerie = QtCharts.QXYSeries
        self._residSerie = QtCharts.QXYSeries

        self._time = {"calculate": 0, "process": 0, "replace": 0}

        self._xMin = self._edDict['experiments'][0][self._first_experiment_name]['_pd_meas_2theta_range_min']
        self._xMax = self._edDict['experiments'][0][self._first_experiment_name]['_pd_meas_2theta_range_max']
        self._xInc = self._edDict['experiments'][0][self._first_experiment_name]['_pd_meas_2theta_range_inc']

        rhochi_calc_chi_sq_by_dictionary(self._cryspyDict, dict_in_out=self._cryspyInOutDict)
        self._xArray = np.arange(self._xMin, self._xMax + self._xInc, self._xInc, dtype=DTYPE)  # -> _xArrayMeasured?
        self._yArrayMeasured = self._cryspyInOutDict[f'pd_{self._first_experiment_name}']['signal_exp'][0]
        self._yArrayMeasuredSigma = self._cryspyInOutDict[f'pd_{self._first_experiment_name}']['signal_exp'][1]
        self._yArray = np.empty_like(self._xArray)  # -> yArrayCalculated?
        self._xArrayProcessed = np.empty_like(self._xArray)
        self._yArrayProcessed = np.empty_like(self._xArray)
        self._yArrayResid = np.empty_like(self._xArray)

        self.wavelengthChanged.connect(self.updateDictData)
        self.resolutionXChanged.connect(self.updateDictData)
        self.stepChanged.connect(self.updateDictData)
        self.radiationProbeChanged.connect(self.updateDictData)
        self.scaleChanged.connect(self.updateDictData)
        self.meas2thetaOffsetChanged.connect(self.updateDictData)

        self.dictDataChanged.connect(self.updateChart)
        self.calculatorChanged.connect(self.updateChart)

    def createEdDict(self, cryspy_obj, cryspy_dict):
        phase_names = [name.replace('crystal_', '') for name in cryspy_dict.keys() if name.startswith('crystal_')]
        experiment_names = [name.replace('pd_', '') for name in cryspy_dict.keys() if name.startswith('pd_')]
        ed_dict = {}
        for data_block in cryspy_obj.items:
            data_block_name = data_block.data_name
            # Phase datablock
            if data_block_name in phase_names:
                ed_dict['phases'] = []
                ed_phase = {data_block_name: {}}
                cryspy_phase = data_block.items
                for item in cryspy_phase:
                    # Space group section
                    if type(item) == cryspy.C_item_loop_classes.cl_2_space_group.SpaceGroup:
                        ed_phase[data_block_name]['_space_group_name_H-M_alt'] = item.name_hm_alt
                    # Cell section
                    elif type(item) == cryspy.C_item_loop_classes.cl_1_cell.Cell:
                        ed_phase[data_block_name]['_cell_length_a'] = item.length_a
                        ed_phase[data_block_name]['_cell_length_b'] = item.length_b
                        ed_phase[data_block_name]['_cell_length_c'] = item.length_c
                        ed_phase[data_block_name]['_cell_angle_alpha'] = item.angle_alpha
                        ed_phase[data_block_name]['_cell_angle_beta'] = item.angle_beta
                        ed_phase[data_block_name]['_cell_angle_gamma'] = item.angle_gamma
                    # Atoms section
                    elif type(item) == cryspy.C_item_loop_classes.cl_1_atom_site.AtomSiteL:
                        ed_atoms = []
                        cryspy_atoms = item.items
                        for cryspy_atom in cryspy_atoms:
                            ed_atom = {}
                            ed_atom['_label'] = cryspy_atom.label
                            ed_atom['_type_symbol'] = cryspy_atom.type_symbol
                            ed_atom['_fract_x'] = cryspy_atom.fract_x
                            ed_atom['_fract_y'] = cryspy_atom.fract_y
                            ed_atom['_fract_z'] = cryspy_atom.fract_z
                            ed_atom['_occupancy'] = cryspy_atom.occupancy
                            ed_atom['_adp_type'] = cryspy_atom.adp_type
                            ed_atom['_B_iso_or_equiv'] = cryspy_atom.b_iso_or_equiv
                            ed_atom['_multiplicity'] = cryspy_atom.multiplicity
                            ed_atom['_Wyckoff_symbol'] = cryspy_atom.wyckoff_symbol
                            ed_atoms.append(ed_atom)
                        ed_phase[data_block_name]['_atom_site'] = ed_atoms
                ed_dict['phases'].append(ed_phase)
                # Experiment datablock
            if data_block_name in experiment_names:
                ed_dict['experiments'] = []
                ed_experiment = {data_block_name: {}}
                cryspy_experiment = data_block.items
                for item in cryspy_experiment:
                    # Ranges section
                    if type(item) == cryspy.C_item_loop_classes.cl_1_range.Range:
                        ed_experiment[data_block_name]['_pd_meas_2theta_range_min'] = item.ttheta_min
                        ed_experiment[data_block_name]['_pd_meas_2theta_range_max'] = item.ttheta_max
                        ed_experiment[data_block_name]['_pd_meas_2theta_range_inc'] = 0.05  # NEED FIX
                    # Setup section
                    elif type(item) == cryspy.C_item_loop_classes.cl_1_setup.Setup:
                        ed_experiment[data_block_name]['_diffrn_radiation_probe'] = item.radiation.replace('neutrons', 'neutron').replace('X-rays', 'x-ray')
                        ed_experiment[data_block_name]['_diffrn_radiation_wavelength'] = item.wavelength
                        ed_experiment[data_block_name]['_pd_meas_2theta_offset'] = item.offset_ttheta
                    # Instrument resolution section
                    elif type(item) == cryspy.C_item_loop_classes.cl_1_pd_instr_resolution.PdInstrResolution:
                        ed_experiment[data_block_name]['_pd_instr_resolution_u'] = item.u
                        ed_experiment[data_block_name]['_pd_instr_resolution_v'] = item.v
                        ed_experiment[data_block_name]['_pd_instr_resolution_w'] = item.w
                        ed_experiment[data_block_name]['_pd_instr_resolution_x'] = item.x
                        ed_experiment[data_block_name]['_pd_instr_resolution_y'] = item.y
                    # Peak assymetries section
                    elif type(item) == cryspy.C_item_loop_classes.cl_1_pd_instr_reflex_asymmetry.PdInstrReflexAsymmetry:
                        ed_experiment[data_block_name]['_pd_instr_reflex_asymmetry_p1'] = item.p1
                        ed_experiment[data_block_name]['_pd_instr_reflex_asymmetry_p2'] = item.p2
                        ed_experiment[data_block_name]['_pd_instr_reflex_asymmetry_p3'] = item.p3
                        ed_experiment[data_block_name]['_pd_instr_reflex_asymmetry_p4'] = item.p4
                    # Phases section
                    elif type(item) == cryspy.C_item_loop_classes.cl_1_phase.PhaseL:
                        ed_phases = []
                        cryspy_phases = item.items
                        for cryspy_phase in cryspy_phases:
                            ed_phase = {}
                            ed_phase['_label'] = cryspy_phase.label
                            ed_phase['_scale'] = cryspy_phase.scale
                            ed_phases.append(ed_phase)
                        ed_experiment[data_block_name]['_phase'] = ed_phases
                    # Background section
                    elif type(item) == cryspy.C_item_loop_classes.cl_1_pd_background.PdBackgroundL:
                        ed_bkg_points = []
                        cryspy_bkg_points = item.items
                        for cryspy_bkg_point in cryspy_bkg_points:
                            ed_bkg_point = {}
                            ed_bkg_point['_2theta'] = cryspy_bkg_point.ttheta
                            ed_bkg_point['_intensity'] = cryspy_bkg_point.intensity
                            ed_bkg_points.append(ed_bkg_point)
                        ed_experiment[data_block_name]['_pd_background'] = ed_bkg_points
                    # Measured data section
                    elif type(item) == cryspy.C_item_loop_classes.cl_1_pd_meas.PdMeasL:
                        ed_meas_points = []
                        cryspy_meas_points = item.items
                        for cryspy_meas_point in cryspy_meas_points:
                            ed_meas_point = {}
                            ed_meas_point['_2theta'] = cryspy_meas_point.ttheta
                            ed_meas_point['_intensity'] = cryspy_meas_point.intensity
                            ed_meas_point['_intensity_sigma'] = cryspy_meas_point.intensity_sigma
                            ed_meas_points.append(ed_meas_point)
                        ed_experiment[data_block_name]['_pd_meas'] = ed_meas_points
                ed_dict['experiments'].append(ed_experiment)
        return ed_dict

    @Property('QVariant', notify=axisRangesChanged)
    def axisRanges(self):
        return self._axisRanges

    @Property(str, notify=calculatorChanged)
    def calculator(self):
        return self._calculator

    @calculator.setter
    def calculator(self, new_value):
        if self._calculator == new_value:
            return
        self._calculator = new_value
        calculateProfile = True
        self.calculatorChanged.emit(calculateProfile)

    @Property(str, notify=radiationProbeChanged)
    def radiationProbe(self):
        return self._radiationProbe

    @radiationProbe.setter
    def radiationProbe(self, new_value):
        if self._radiationProbe == new_value:
            return
        self._radiationProbe = new_value
        calculateProfile = True
        self.radiationProbeChanged.emit(calculateProfile)

    @Property(int, notify=arraySizeChanged)
    def arraySize(self):
        return self._arraySize

    @arraySize.setter
    def arraySize(self, new_value):
        if self._arraySize == new_value:
            return
        self._arraySize = new_value
        self.arraySizeChanged.emit()

    @Property(float, notify=wavelengthChanged)
    def wavelength(self):
        return self._wavelength

    @wavelength.setter
    def wavelength(self, new_value):
        if self._wavelength == new_value:
            return
        self._wavelength = new_value
        calculateProfile = True
        self.wavelengthChanged.emit(calculateProfile)

    @Property(float, notify=scaleChanged)
    def scale(self):
        return self._scale

    @scale.setter
    def scale(self, new_value):
        if self._scale == new_value:
            return
        self._scale = new_value
        calculateProfile = False
        self.scaleChanged.emit(calculateProfile)

    @Property(float, notify=meas2thetaOffsetChanged)
    def meas2thetaOffset(self):
        return self._meas2thetaOffset

    @meas2thetaOffset.setter
    def meas2thetaOffset(self, new_value):
        if self._meas2thetaOffset == new_value:
            return
        self._meas2thetaOffset = new_value
        calculateProfile = False
        self.meas2thetaOffsetChanged.emit(calculateProfile)

    @Property(float, notify=stepChanged)
    def step(self):
        return self._step

    @step.setter
    def step(self, new_value):
        if self._step == new_value:
            return
        self._step = new_value
        calculateProfile = True
        self.stepChanged.emit(calculateProfile)

    @Property(float, notify=resolutionXChanged)
    def resolutionX(self):
        return self._resolutionX

    @resolutionX.setter
    def resolutionX(self, new_value):
        if self._resolutionX == new_value:
            return
        self._resolutionX = new_value
        calculateProfile = True
        self.resolutionXChanged.emit(calculateProfile)

    def updateDictData(self, calculateProfile):
        self._edDict['experiments'][0][self._first_experiment_name]['_diffrn_radiation_wavelength'] = self._wavelength
        self._edDict['experiments'][0][self._first_experiment_name]['_pd_instr_resolution_x'] = self._resolutionX
        self._edDict['experiments'][0][self._first_experiment_name]['_pd_meas_2theta_range_inc'] = self._step
        self._edDict['experiments'][0][self._first_experiment_name]['_diffrn_radiation_probe'] = self._radiationProbe

        self._edDict['experiments'][0][self._first_experiment_name]['_pd_meas_2theta_offset'] = self._meas2thetaOffset
        self._edDict['experiments'][0][self._first_experiment_name]['_phase'][0]['_scale'] = self._scale

        if self._calculator == 'crysfml':
            self._crysfmlDict['experiments'][0][self._first_experiment_name]['_diffrn_radiation_wavelength'] = self._wavelength
            self._crysfmlDict['experiments'][0][self._first_experiment_name]['_pd_instr_resolution_x'] = self._resolutionX
            self._crysfmlDict['experiments'][0][self._first_experiment_name]['_pd_meas_2theta_range_inc'] = self._step
            self._crysfmlDict['experiments'][0][self._first_experiment_name]['_diffrn_radiation_probe'] = self._radiationProbe

        elif self._calculator == 'cryspy':
            self._cryspyDict[f'pd_{self._first_experiment_name}']['wavelength'][0] = self._wavelength
            self._cryspyDict[f'pd_{self._first_experiment_name}']['resolution_parameters'][3] = self._resolutionX
            self._cryspyDict[f'pd_{self._first_experiment_name}']['radiation'][0] = self._radiationProbe.replace('neutron', 'neutrons').replace('x-ray', 'X-rays')

        self.dictDataChanged.emit(calculateProfile)

    @Property('QVariant', notify=timeChanged)
    def time(self):
        return self._time

    @time.setter
    def time(self, new_value):
        if self._time == new_value:
            return
        self._time = new_value
        self.timeChanged.emit()

    @Property('QVariant', notify=measSerieChanged)
    def measSerie(self):
        return self._measSerie

    @measSerie.setter
    def measSerie(self, new_serie):
        if self._measSerie == new_serie:
            return
        self._measSerie = new_serie

        self.replaceChartMeasData()
        self.updateChartAxisRanges()

        self.measSerieChanged.emit()

    @Property('QVariant', notify=calcSerieChanged)
    def calcSerie(self):
        return self._calcSerie

    @calcSerie.setter
    def calcSerie(self, new_serie):
        if self._calcSerie == new_serie:
            return
        self._calcSerie = new_serie
        self.calcSerieChanged.emit()

    #
    @Property('QVariant', notify=residSerieChanged)
    def residSerie(self):
        return self._residSerie

    @residSerie.setter
    def residSerie(self, new_serie):
        if self._residSerie == new_serie:
            return
        self._residSerie = new_serie
        self.residSerieChanged.emit()

    def timing(self, f):
        return timeit.timeit(stmt=f, number=1)

    def calculateDiffractionPattern(self):
        if self._calculator == 'crysfml':
            self._yArray = CFML_api.powder_mod.simulation(self._crysfmlDict)[0].astype(np.float64)
        elif self._calculator == 'cryspy':
            rhochi_calc_chi_sq_by_dictionary(self._cryspyDict,
                                             dict_in_out=self._cryspyInOutDict,
                                             flag_use_precalculated_data=False,
                                             flag_calc_analytical_derivatives=False)
            self._yArray = self._cryspyInOutDict[f'pd_{self._first_experiment_name}']['signal_plus']
        self.arraySize = int(self._xArray.size)

    def processDiffractionPattern(self):
        # params
        scale = DTYPE(self._scale)
        offset = DTYPE(self._meas2thetaOffset)
        # processed x-array
        self._xArrayProcessed = self._xArray + offset
        # background y-array
        background = self._edDict['experiments'][0][self._first_experiment_name]['_pd_background']
        xp = np.array([item['_2theta'] for item in background], dtype=DTYPE)
        fp = np.array([item['_intensity'] for item in background], dtype=DTYPE)
        backgroundYArray = np.interp(self._xArrayProcessed, xp, fp)
        # processed y-array
        multiplier = scale * DTYPE(100.0) / self._yArray.max()
        self._yArrayProcessed = self._yArray * multiplier + backgroundYArray
        # residual y-array
        self._yArrayResid = self._yArrayMeasured - self._yArrayProcessed

    def replaceChartMeasData(self):
        print('--------------- self._measSerie', self._measSerie)
        self._measSerie.replaceNp(self._xArray, self._yArrayMeasured)
        self.measSerieChanged.emit()

    def replaceChartCalcData(self):
        self._calcSerie.replaceNp(self._xArray, self._yArrayProcessed)
        self._residSerie.replaceNp(self._xArray, self._yArrayResid)
        self.calcSerieChanged.emit()
        self.residSerieChanged.emit()

    def updateChartAxisRanges(self):
        self._axisRanges = {'xMin': float(self._xArray.min()),
                            'xMax': float(self._xArray.max()),
                            'yMin': float(self._yArrayMeasured.min()),
                            'yMax': float(self._yArrayMeasured.max())}
        self.axisRangesChanged.emit()

    def chisq_array(self):
        return np.square((self._yArrayProcessed - self._yArrayMeasured) / self._yArrayMeasuredSigma)

    def chisq_sum(self):
        return np.sum(self.chisq_array())

    # FIT 1: scipy optimize

    def processFuncScipy(self, params):
        self._scale = params[0]
        print(f'scale: {self._scale}')
        self.processDiffractionPattern()
        return self.chisq_sum()

    @Slot()
    def fitProcessedOnlyScipy(self):
        fun = self.processFuncScipy
        x0 = np.array([self._scale])
        method = 'lm' # ‘trf’, ‘dogbox’, ‘lm’
        result = sc.optimize.least_squares(fun, x0, method=method)
        print(result)
        if result.success:
            self._scale = result.x[0]
            calculateProfile = False
            self.scaleChanged.emit(calculateProfile)

    def calculateAndProcessFuncScipy(self, params, xArray, yArray):
        self._wavelength = params[0]
        self._resolutionX = params[1]
        print(f'wavelength: {self._wavelength}, resolutionX: {self._resolutionX}')
        if self._calculator == 'crysfml':
            self._crysfmlDict['experiments'][0][self._first_experiment_name]['_diffrn_radiation_wavelength'] = self._wavelength
            self._crysfmlDict['experiments'][0][self._first_experiment_name]['_pd_instr_resolution_x'] = self._resolutionX
        elif self._calculator == 'cryspy':
            self._cryspyDict[f'pd_{self._first_experiment_name}']['wavelength'][0] = self._wavelength
            self._cryspyDict[f'pd_{self._first_experiment_name}']['resolution_parameters'][3] = self._resolutionX
        self.calculateDiffractionPattern()
        self.processDiffractionPattern()
        return self.chisq_array()

    @Slot()
    def fitWithCalcScipy(self):
        fun = self.calculateAndProcessFuncScipy
        x0 = np.array([self._wavelength, self._resolutionX])
        method = 'lm' # ‘trf’, ‘dogbox’, ‘lm’
        result = sc.optimize.least_squares(fun, x0, args=(self._xArray, self._yArray), method=method)
        print(result)
        if result.success:
            self._wavelength = result.x[0]
            self._resolutionX = result.x[1]
            calculateProfile = True
            self.wavelengthChanged.emit(calculateProfile)
            self.resolutionXChanged.emit(calculateProfile)

    # FIT 1: lmfit minimize

    def processFuncLmfit(self, params):
        self._scale = params['scale'].value
        print(f'scale: {self._scale}')
        self.processDiffractionPattern()
        return self.chisq_array()

    @Slot()
    def fitProcessedOnlyLmfit(self):
        params = lmfit.Parameters()
        #params.add('scale', value=self._scale, min=1.000001, max=1000.0)
        params.add('scale', value=self._scale)
        result = lmfit.minimize(self.processFuncLmfit, params, method='least_squares')
        lmfit.report_fit(result)
        if result.success:
            self._scale = result.params['scale'].value
            calculateProfile = False
            self.scaleChanged.emit(calculateProfile)

    def calculateAndProcessFuncLmfit(self, params):
        self._wavelength = params['wavelength'].value
        self._resolutionX = params['resolutionX'].value
        print(f'wavelength: {self._wavelength}, resolutionX: {self._resolutionX}')
        if self._calculator == 'crysfml':
            self._crysfmlDict['experiments'][0][self._first_experiment_name]['_diffrn_radiation_wavelength'] = self._wavelength
            self._crysfmlDict['experiments'][0][self._first_experiment_name]['_pd_instr_resolution_x'] = self._resolutionX
        elif self._calculator == 'cryspy':
            self._cryspyDict[f'pd_{self._first_experiment_name}']['wavelength'][0] = self._wavelength
            self._cryspyDict[f'pd_{self._first_experiment_name}']['resolution_parameters'][3] = self._resolutionX
        self.calculateDiffractionPattern()
        self.processDiffractionPattern()
        return self.chisq_array()

    @Slot()
    def fitWithCalcLmfit(self):
        params = lmfit.Parameters()
        #params.add('wavelength', value=self._wavelength, min=1.89, max=1.93)
        params.add('wavelength', value=self._wavelength)
        params.add('resolutionX', value=self._resolutionX)
        result = lmfit.minimize(self.calculateAndProcessFuncLmfit, params, method='least_squares')
        lmfit.report_fit(result)
        if result.success:
            self._wavelength = result.params['wavelength'].value
            self._resolutionX = result.params['resolutionX'].value
            calculateProfile = True
            self.wavelengthChanged.emit(calculateProfile)
            self.resolutionXChanged.emit(calculateProfile)

    # FIT: Minuit

    def processFuncMinuit(self, scale):
        self._scale = scale
        print(f'scale: {self._scale}')
        self.processDiffractionPattern()
        return self.chisq_sum()

    @Slot()
    def fitProcessedOnlyMinuit(self):
        m = iminuit.Minuit(self.processFuncMinuit, self._scale)
        #m.limits = [(0.000001, 1000.0)]
        #m.errors = (100.0)
        m.errordef = iminuit.Minuit.LEAST_SQUARES
        #m.simplex()
        m.migrad()  # finds minimum of least_squares function
        m.hesse()   # accurately computes uncertainties
        print(m.params)
        for p, v, e in zip(m.parameters, m.values, m.errors):
            print(f"{p} = {v} +- {e}")
        if True:
            self._scale = m.values[0]
            calculateProfile = False
            self.scaleChanged.emit(calculateProfile)

    def calculateAndProcessFuncMinuit(self, wavelength, resolutionX):
        self._wavelength = wavelength
        self._resolutionX = resolutionX
        print(f'wavelength: {self._wavelength}, resolutionX: {self._resolutionX}')
        if self._calculator == 'crysfml':
            self._crysfmlDict['experiments'][0][self._first_experiment_name]['_diffrn_radiation_wavelength'] = self._wavelength
            self._crysfmlDict['experiments'][0][self._first_experiment_name]['_pd_instr_resolution_x'] = self._resolutionX
        elif self._calculator == 'cryspy':
            self._cryspyDict[f'pd_{self._first_experiment_name}']['wavelength'][0] = self._wavelength
            self._cryspyDict[f'pd_{self._first_experiment_name}']['resolution_parameters'][3] = self._resolutionX
        self.calculateDiffractionPattern()
        self.processDiffractionPattern()
        return self.chisq_sum()

    @Slot()
    def fitWithCalcMinuit(self):
        m = iminuit.Minuit(self.calculateAndProcessFuncMinuit, self._wavelength, self._resolutionX)
        #m.limits = [(1.89, 1.93)]
        #m.errors = (0.2)
        m.errordef = iminuit.Minuit.LEAST_SQUARES
        #m.simplex()
        m.migrad()  # finds minimum of least_squares function
        m.hesse()   # accurately computes uncertainties
        print(m.params)
        for p, v, e in zip(m.parameters, m.values, m.errors):
            print(f"{p} = {v} +- {e}")
        if True:
            self._wavelength = m.values[0]
            self._resolutionX = m.values[1]
            calculateProfile = True
            self.wavelengthChanged.emit(calculateProfile)
            self.resolutionXChanged.emit(calculateProfile)

    ##################

    @Slot(bool)
    def updateChart(self, calculateProfile):
        if calculateProfile:
            self._time['calculate'] = self.timing(self.calculateDiffractionPattern)
        else:
            self._time['calculate'] = 0
        self._time['process'] = self.timing(self.processDiffractionPattern)
        self._time['replace'] = self.timing(self.replaceChartCalcData)
        self.timeChanged.emit()


if __name__ == '__main__':
    os.environ['QSG_RHI_BACKEND'] = 'opengl'  # For QtCharts LineSeries useOpenGL
    os.environ['FULLPROF'] = os.path.abspath(os.path.dirname(__file__))
    app = QApplication(sys.argv)
    QSurfaceFormat.setDefaultFormat(QQuick3D.idealSurfaceFormat(4))  # QtQuick3D
    proxy = Proxy()
    engine = QQmlApplicationEngine()
    engine.rootContext().setContextProperty("proxy", proxy)
    engine.load("main.qml")
    if not engine.rootObjects():
        sys.exit(-1)
    sys.exit(app.exec())
