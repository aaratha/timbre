import QtQuick 6.4
import QtQuick.Controls 6.4

ApplicationWindow {
    width: 640
    height: 480
    visible: true
    color: "black"

    property int rectCount: Math.max(1, ui ? ui.binCount : 1)
    property int minRectWidth: 8
    property int activeIndex: -1
    property int pointSize: 6
    property int plotPadding: 16

    Item {
        id: binStrip
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: parent.top
        anchors.margins: 16
        height: 80

        Row {
            id: row
            anchors.fill: parent
            spacing: 4

            property real cellWidth: Math.max(minRectWidth, (row.width - (rectCount - 1) * row.spacing) / rectCount)

            Repeater {
                model: rectCount
                delegate: Rectangle {
                    width: row.cellWidth
                    height: row.height
                    color: index === activeIndex ? "#f39c12" : "#e74c3c"
                    radius: 6
                }
            }
        }

        MouseArea {
            anchors.fill: parent
            hoverEnabled: false

            function updateIndex(xPos) {
                var cell = row.cellWidth + row.spacing
                var idx = Math.floor(xPos / cell)
                if (idx < 0 || idx >= rectCount) {
                    return
                }
                if (activeIndex !== idx) {
                    activeIndex = idx
                    ui.rectangleClicked(idx)
                }
            }

            onPressed: function(mouse) { updateIndex(mouse.x) }
            onPositionChanged: function(mouse) { if (pressed) updateIndex(mouse.x) }
        }
    }

    Rectangle {
        id: plotArea
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: binStrip.bottom
        anchors.bottom: parent.bottom
        anchors.margins: 16
        color: "#101418"
        border.color: "#2c3e50"
        radius: 8

        Repeater {
            model: ui ? ui.umapPoints : []
            delegate: Item {
                width: pointSize
                height: pointSize

                readonly property real xNorm: modelData.x
                readonly property real yNorm: modelData.y

                x: plotPadding + xNorm * (plotArea.width - 2 * plotPadding) - width / 2
                y: plotPadding + (1 - yNorm) * (plotArea.height - 2 * plotPadding) - height / 2

                Rectangle {
                    anchors.fill: parent
                    radius: width / 2
                    color: modelData.index === activeIndex ? "#f39c12" : "#1abc9c"
                }
            }
        }

        MouseArea {
            anchors.fill: parent
            hoverEnabled: false

            function updateIndex(xPos, yPos) {
                if (!ui || !ui.umapPoints || ui.umapPoints.length === 0) {
                    return
                }
                var bestIndex = -1
                var bestDist = 1e9
                var w = plotArea.width - 2 * plotPadding
                var h = plotArea.height - 2 * plotPadding
                for (var i = 0; i < ui.umapPoints.length; ++i) {
                    var p = ui.umapPoints[i]
                    var px = plotPadding + p.x * w
                    var py = plotPadding + (1 - p.y) * h
                    var dx = xPos - px
                    var dy = yPos - py
                    var d2 = dx * dx + dy * dy
                    if (d2 < bestDist) {
                        bestDist = d2
                        bestIndex = p.index
                    }
                }
                if (bestIndex >= 0 && activeIndex !== bestIndex) {
                    activeIndex = bestIndex
                    ui.rectangleClicked(bestIndex)
                }
            }

            onPressed: function(mouse) { updateIndex(mouse.x, mouse.y) }
            onPositionChanged: function(mouse) { if (pressed) updateIndex(mouse.x, mouse.y) }
        }
    }
}
