import QtQuick
import QtQuick.Controls
import QtQuick.Controls.Universal
import QtQuick3D
import QtCharts

ApplicationWindow {
    id: applicationWindow

    property int margins: 10
    property int labelHeight: 24
    property int editWidth: 70
    property int sliderWidth: 200
    property bool useOpenGL: true

    Universal.theme: Universal.Light
    Universal.accent: Universal.Cyan

    visible: true
    width: minimumWidth
    height: minimumHeight
    minimumWidth: 1600
    minimumHeight: 750

    Row {
        id: layout

        Column {

            // Chart

            ChartView {
                id: chartView
                width: applicationWindow.width - sidebar.width
                height: applicationWindow.height * 0.7
                antialiasing: true
                legend.visible: false
                title: 'Calculated diffraction pattern PbSO4'
                titleColor: '#999'
                //animationOptions: ChartView.SeriesAnimations
                //animationDuration: 100
                ScatterSeries {
                    id: measSerie
                    axisX: axisX
                    axisY: axisY
                    useOpenGL: false
                    markerSize: 6
                    borderWidth: 1
                    color: 'cornflowerblue'
                    borderColor: measSerie.color
                }
                LineSeries {
                    id: calcSerie
                    axisX: axisX
                    axisY: axisY
                    useOpenGL: false
                    color: 'coral'
                    width: 2
                }
                ValueAxis {
                    id: axisX
                    min: proxy.axisRanges.xMin
                    max: proxy.axisRanges.xMax
                }
                ValueAxis {
                    id: axisY
                    min: proxy.axisRanges.yMin
                    max: proxy.axisRanges.yMax
                }
                Component.onCompleted: {

                    proxy.measSerie = measSerie
                    print('+++++++++++++++', measSerie, proxy.measSerie)
                    proxy.calcSerie = calcSerie
                    const calculateProfile = true
                    proxy.updateChart(calculateProfile)
                }
                // Zoom rectangle
                Rectangle{
                    id: recZoom
                    property int xScaleZoom: 0
                    property int yScaleZoom: 0
                    visible: false
                    transform: Scale { origin.x: 0; origin.y: 0; xScale: recZoom.xScaleZoom; yScale: recZoom.yScaleZoom}
                    border.color: "#999"
                    border.width: 1
                    opacity: 0.9
                    color: "transparent"
                    Rectangle {
                        anchors.fill: parent
                        opacity: 0.5
                        color: recZoom.border.color
                    }
                }
                // Left mouse button events
                MouseArea {
                    anchors.fill: chartView
                    acceptedButtons: Qt.LeftButton
                    onPressed: {
                        recZoom.x = mouseX
                        recZoom.y = mouseY
                        recZoom.visible = true
                    }
                    onMouseXChanged: {
                        if (mouseX > recZoom.x) {
                            recZoom.xScaleZoom = 1
                            recZoom.width = Math.min(mouseX, chartView.width) - recZoom.x
                        } else {
                            recZoom.xScaleZoom = -1
                            recZoom.width = recZoom.x - Math.max(mouseX, 0)
                        }
                    }
                    onMouseYChanged: {
                        if (mouseY > recZoom.y) {
                            recZoom.yScaleZoom = 1
                            recZoom.height = Math.min(mouseY, chartView.height) - recZoom.y
                        } else {
                            recZoom.yScaleZoom = -1
                            recZoom.height = recZoom.y - Math.max(mouseY, 0)
                        }
                    }
                    onReleased: {
                        const x = Math.min(recZoom.x, mouseX) - chartView.anchors.leftMargin
                        const y = Math.min(recZoom.y, mouseY) - chartView.anchors.topMargin
                        const width = recZoom.width
                        const height = recZoom.height
                        chartView.zoomIn(Qt.rect(x, y, width, height))
                        recZoom.visible = false
                    }
                }
                // Right mouse button events
                MouseArea {
                    anchors.fill: chartView
                    acceptedButtons: Qt.RightButton
                    onClicked: chartView.zoomReset()
                }
            }


            ChartView {
                id: chartViewResid
                width: chartView.width
                height: applicationWindow.height - chartView.height
                antialiasing: true
                legend.visible: false
                LineSeries {
                    id: residSerie
                    axisX: axisXResid
                    axisY: axisYResid
                    useOpenGL: false
                    color: 'gray'
                    width: 2
                }
                ValueAxis {
                    id: axisXResid
                    min: axisX.min
                    max: axisX.max
                }
                ValueAxis {
                    id: axisYResid
                    min: -1000 //proxy.axisRanges.yMin
                    max: 1000 //proxy.axisRanges.yMax
                }
                Component.onCompleted: {
                    proxy.residSerie = residSerie
                }


            }


        }



        // Parameters

        Column {
            id: sidebar

            padding: margins
            topPadding: margins * 3.5
            leftPadding: 0
            spacing: margins

            Label {
                text: 'Parameters'
                color: '#999'
            }

            // calculator
            Row {
                Label {
                    height: labelHeight * 1.3
                    verticalAlignment: Text.AlignVCenter
                    text: 'Calculator    '
                }
                ButtonGroup {
                    id: calculatorRadioGroup
                    onClicked: (button) => proxy.calculator = button.text
                }
                RadioButton {
                    checked: proxy.calculator === text
                    text: 'crysfml'
                    ButtonGroup.group: calculatorRadioGroup
                }
                RadioButton {
                    checked: proxy.calculator === text
                    text: 'cryspy'
                    ButtonGroup.group: calculatorRadioGroup
                }
            }

            // _diffrn_radiation_probe
            Row {
                Label {
                    height: labelHeight * 1.3
                    verticalAlignment: Text.AlignVCenter
                    text: 'Radiation    '
                }
                ButtonGroup {
                    id: probeRadioGroup
                    onClicked: (button) => proxy.radiationProbe = button.text
                }
                RadioButton {
                    checked: proxy.radiationProbe === text
                    text: 'neutron'
                    ButtonGroup.group: probeRadioGroup
                }
                RadioButton {
                    checked: proxy.radiationProbe === text
                    text: 'x-ray'
                    ButtonGroup.group: probeRadioGroup
                }
            }

            // _diffrn_radiation_wavelength
            Row {
                spacing: margins
                Label {
                    height: labelHeight
                    verticalAlignment: Text.AlignVCenter
                    text: 'wavelength'
                }
                TextField {
                    width: editWidth
                    text: proxy.wavelength.toFixed(4)
                    onEditingFinished: proxy.wavelength = parseFloat(text)
                }
                Slider {
                    width: sliderWidth
                    from: 1.89
                    to: 1.93
                    value: proxy.wavelength
                    onMoved: {
                        if (useOpenGL) calcSerie.useOpenGL = true
                        proxy.wavelength = value
                        if (useOpenGL) disableOpenGLTimer.restart()
                    }
                }
            }

            // _pd_instr_resolution_y
            Row {
                spacing: margins
                Label {
                    height: labelHeight
                    verticalAlignment: Text.AlignVCenter
                    text: 'resolution x'
                }
                TextField {
                    width: editWidth
                    text: proxy.resolutionX.toFixed(3)
                    onEditingFinished: proxy.resolutionX = parseFloat(text)
                }
                Slider {
                    width: sliderWidth
                    from: 0.0
                    to: 0.5
                    value: proxy.resolutionX
                    onMoved: {
                        if (useOpenGL) calcSerie.useOpenGL = true
                        proxy.resolutionX = value
                        if (useOpenGL) disableOpenGLTimer.restart()
                    }
                }
            }

            // _pd_meas_2theta_offset
            Row {
                enabled: false
                spacing: margins
                Label {
                    height: labelHeight
                    verticalAlignment: Text.AlignVCenter
                    text: '2theta offset'
                }
                TextField {
                    width: editWidth
                    text: proxy.meas2thetaOffset.toFixed(3)
                    onEditingFinished: proxy.meas2thetaOffset = parseFloat(text)
                }
                Slider {
                    width: sliderWidth
                    from: -0.5
                    to: 0.5
                    value: proxy.meas2thetaOffset
                    onMoved: {
                        if (useOpenGL) calcSerie.useOpenGL = true
                        proxy.meas2thetaOffset = value
                        if (useOpenGL) disableOpenGLTimer.restart()
                    }
                }
            }

            // _scale
            Row {
                spacing: margins
                Label {
                    height: labelHeight
                    verticalAlignment: Text.AlignVCenter
                    text: 'phase scale'
                }
                TextField {
                    width: editWidth
                    text: proxy.scale.toFixed(2)
                    onEditingFinished: proxy.scale = parseFloat(text)
                }
                Slider {
                    width: sliderWidth
                    from: 1
                    to: 100
                    value: proxy.scale
                    onMoved: {
                        if (useOpenGL) calcSerie.useOpenGL = true
                        proxy.scale = value
                        if (useOpenGL) disableOpenGLTimer.restart()
                    }
                }
            }

            // _pd_meas_2theta_range_inc
            Row {
                enabled: proxy.calculator === 'crysfml'
                spacing: margins
                Label {
                    height: labelHeight
                    verticalAlignment: Text.AlignVCenter
                    text: '2theta step'
                }
                TextField {
                    width: editWidth
                    text: proxy.step.toFixed(4)
                    onEditingFinished: proxy.step = parseFloat(text)
                }
                Slider {
                    width: sliderWidth
                    from: 0.001
                    to: 1.0
                    value: proxy.step
                    onMoved: {
                        if (useOpenGL) calcSerie.useOpenGL = true
                        proxy.step = value
                        if (useOpenGL) disableOpenGLTimer.restart()
                    }
                }
            }

            // Measured data
            CheckBox {
                text: 'Hide measured data'
                checked: measSerie.visible
                onCheckedChanged: measSerie.visible = checked
            }

            // OpenGL
            CheckBox {
                text: 'Use OpenGL during chart update'
                checked: useOpenGL
                onCheckedChanged: useOpenGL = checked
            }

            // Status

            Column {
                spacing: margins
                Label {
                    text: 'Status'
                    color: '#999'
                }
                Label {
                    text: `Data array size: ${proxy.arraySize} points`
                }
                Label {
                    text: `Pattern calculation time: ${proxy.time.calculate.toFixed(4)} sec`
                }
                Label {
                    text: `Pattern post-processing time: ${proxy.time.process.toFixed(4)} sec`
                }
                Label {
                    text: `Chart data replacement time: ${proxy.time.replace.toFixed(4)} sec`
                }
            }

            // Fit

            Row {
                spacing: margins
                Button {
                    text: 'Fit scale: scipy.optimize.least_squares lm'
                    onClicked: proxy.fitProcessedOnlyScipy()
                }
                Button {
                    text: 'Fit wl.+res.x: scipy.optimize.least_squares lm'
                    onClicked: proxy.fitWithCalcScipy()
                }
            }

            Row {
                spacing: margins
                Button {
                    text: 'Fit scale: lmfit.minimize least_squares'
                    onClicked: proxy.fitProcessedOnlyLmfit()
                }
                Button {
                    text: 'Fit wl.+res.x: lmfit.minimize least_squares'
                    onClicked: proxy.fitWithCalcLmfit()
                }
            }

            Row {
                spacing: margins
                Button {
                    text: 'Fit scale: minuit migrad+hesse'
                    onClicked: proxy.fitProcessedOnlyMinuit()
                }
                Button {
                    text: 'Fit wl.+res.x: minuit migrad+hesse'
                    onClicked: proxy.fitWithCalcMinuit()
                }
            }


        }
    }

    // Use OpenGL on slider move only
    Timer {
        id: disableOpenGLTimer
        interval: 500
        onTriggered: calcSerie.useOpenGL = false
    }

}
