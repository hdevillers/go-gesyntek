package heatmap

import (
	"bufio"
	"errors"
	"fmt"
	"image/color"
	"math"
	"os"
	"strconv"
	"strings"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/palette"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
	"gonum.org/v1/plot/vg/vgpdf"
)

type HeatMapGrid struct {
	grid       [][]float64
	N          int
	M          int
	resolution float64
	minX       float64
	minY       float64
}

func (hmp HeatMapGrid) Dims() (c, r int) {
	return hmp.N, hmp.M
}

func (hmp HeatMapGrid) X(c int) float64 {
	return hmp.minX + float64(c)*hmp.resolution
}

func (hmp HeatMapGrid) Y(r int) float64 {
	return hmp.minY + float64(r)*hmp.resolution
}

func (hmp HeatMapGrid) Z(c, r int) float64 {
	return hmp.grid[c][r]
}

/*
CLASS: HeatMapTicks
*/
type HeatMapTicks struct {
	Labels     []string
	Resolution float64
	Offset     string
}

func NewHeatMapTicks(l *[]string, r float64, o int) *HeatMapTicks {
	var hmt HeatMapTicks
	hmt.Labels = *l
	hmt.Resolution = r
	hmt.Offset = ""
	for range o {
		hmt.Offset += " "
	}
	return &hmt
}

func (hmt HeatMapTicks) Ticks(min, max float64) []plot.Tick {
	var pt []plot.Tick
	/*for i := math.Trunc(min); i <= max; i++ {
		at := min + float64(i)*hmt.Resolution
		pt = append(pt, plot.Tick{Value: at, Label: hmt.Labels[int(i)]})
	}*/
	for i := range len(hmt.Labels) {
		at := float64(i) * hmt.Resolution
		pt = append(pt, plot.Tick{Value: at, Label: hmt.Labels[int(i)] + hmt.Offset})
	}
	return pt
}

/*
CLASS: HeatMapDistance
*/
type HeatMapDistance struct {
	Data       [][]float64
	Labels     []string
	XLabOffset int
	YLabOffset int
	Dim        int
	UpMax      float64
	DoMax      float64
	NDiv       int
	Margin     float64
	Title      string
	TitleSize  float64
	Width      float64
	Height     float64
	Output     string
}

func NewHeatMapDistance(n, xo, yo int, m float64, t string, ts float64, w, h float64, o string) *HeatMapDistance {
	var hmd HeatMapDistance
	hmd.UpMax = 0.0
	hmd.DoMax = 0.0
	hmd.Dim = 0
	hmd.NDiv = n
	hmd.XLabOffset = xo
	hmd.YLabOffset = yo
	hmd.Margin = m
	hmd.Title = t
	hmd.TitleSize = ts
	hmd.Width = w
	hmd.Height = h
	hmd.Output = o
	return &hmd
}

func (hmd *HeatMapDistance) LoadDistanceFile(fd string) error {
	// Default value
	noVal := math.NaN()

	// Create file handler
	fh, err := os.Open(fd)
	if err != nil {
		return err
	}
	defer fh.Close()

	// Create the buffer
	fb := bufio.NewScanner(fh)

	// Create temporary variables to contain data
	upDist := make(map[string]float64, 0)
	doDist := make(map[string]float64, 0)
	upOk := 0
	doOk := 0
	keys := make(map[string]bool, 0)
	hmd.Labels = make([]string, 0)

	// Skip the header
	fb.Scan()
	err = fb.Err()
	if err != nil {
		return err
	}

	// Scan each line
	for fb.Scan() {
		err = fb.Err()
		if err != nil {
			return err
		}

		// Read line
		line := fb.Text()
		elem := strings.Split(line, "\t")

		// Store new keys if necessary
		_, ok := keys[elem[0]]
		if !ok {
			keys[elem[0]] = true
			hmd.Labels = append(hmd.Labels, elem[0])
		}
		_, ok = keys[elem[1]]
		if !ok {
			keys[elem[1]] = true
			hmd.Labels = append(hmd.Labels, elem[1])
		}

		// Store values
		cKey := elem[0] + elem[1]
		if elem[2] == "NA" {
			upDist[cKey] = noVal
		} else {
			upDist[cKey], err = strconv.ParseFloat(elem[2], 64)
			if err != nil {
				return err
			}
			upOk++
		}
		if elem[3] == "NA" {
			doDist[cKey] = noVal
		} else {
			doDist[cKey], err = strconv.ParseFloat(elem[3], 64)
			if err != nil {
				return err
			}
			doOk++
		}
	}

	// Check
	if upOk == 0 {
		return errors.New("no distance values for upstream comparison")
	}
	if doOk == 0 {
		return errors.New("no distance values for downstream comparison")
	}

	// Prepare the square matrix
	nKeys := len(hmd.Labels)
	hmd.Dim = nKeys
	hmd.Data = make([][]float64, nKeys)
	for i := range nKeys {
		hmd.Data[i] = make([]float64, nKeys)
	}

	// Initialize min/max
	hmd.UpMax = float64(0.0)
	hmd.DoMax = float64(0.0)
	for i, j := 0, true; j && i < nKeys-1; {
		key := hmd.Labels[i] + hmd.Labels[i+1]
		val := upDist[key]
		if !math.IsNaN(val) {
			//upMin = val
			hmd.UpMax = val
			j = false
		}
	}
	for i, j := 0, true; j && i < nKeys-1; {
		key := hmd.Labels[i] + hmd.Labels[i+1]
		val := doDist[key]
		if !math.IsNaN(val) {
			//doMin = val
			hmd.DoMax = val
			j = false
		}
	}

	// Fill the square matrix
	for i := range nKeys - 1 {
		for j := i + 1; j < nKeys; j++ {
			key := hmd.Labels[i] + hmd.Labels[j]
			upVal := upDist[key]
			doVal := doDist[key]
			hmd.Data[i][j] = upVal
			hmd.Data[j][i] = doVal
			if !math.IsNaN(upVal) && upVal > hmd.UpMax {
				hmd.UpMax = upVal
			}
			if !math.IsNaN(doVal) && doVal > hmd.DoMax {
				hmd.DoMax = doVal
			}
		}
	}
	return nil
}

func (hmd *HeatMapDistance) Plot() error {
	if hmd.Dim == 0 {
		return errors.New("cannot plot distance heatmap: no data loaded")
	}
	// Create the grid
	grid := HeatMapGrid{
		grid:       hmd.Data,
		N:          hmd.Dim,
		M:          hmd.Dim,
		minX:       0.0,
		minY:       0.0,
		resolution: 0.5,
	}

	// Create the palette
	pal := palette.Heat(hmd.NDiv, 1.0)
	hea := plotter.NewHeatMap(grid, pal)
	hea.Min = 0.0
	hea.Max = math.Max(hmd.UpMax, hmd.DoMax) //* 1.1 // Max value is ignored
	hea.NaN = color.RGBA{150, 150, 150, 255}

	// New plot, add the heatmap
	plo := plot.New()
	plo.Add(hea)

	// Axes
	xhmt := NewHeatMapTicks(&hmd.Labels, grid.resolution, hmd.XLabOffset)
	yhmt := NewHeatMapTicks(&hmd.Labels, grid.resolution, hmd.YLabOffset)
	plo.X.Tick.Marker = xhmt
	plo.Y.Tick.Marker = yhmt
	plo.X.Tick.Label.Rotation = math.Pi / 2
	plo.X.Tick.Label.XAlign = draw.XRight
	plo.X.Tick.Label.YAlign = draw.YCenter

	// Legend
	leg := plot.NewLegend()
	thumbs := plotter.PaletteThumbnailers(pal)
	for i := len(thumbs) - 1; i >= 0; i-- {
		t := thumbs[i]
		if i != 0 && i != len(thumbs)-1 {
			leg.Add("", t)
			continue
		}
		var val float64
		switch i {
		case 0:
			val = hea.Min
		case len(thumbs) - 1:
			val = hea.Max
		}
		leg.Add(fmt.Sprintf("%.2g", val), t)
	}
	leg.Top = false

	// Create the canvas linked to the pdf
	img := vgpdf.New(vg.Length(hmd.Width)*vg.Inch, vg.Length(hmd.Height)*vg.Inch)
	drc := draw.New(img)
	drc = draw.Crop(drc,
		vg.Length(hmd.Margin)*vg.Inch,
		-vg.Length(hmd.Margin)*vg.Inch,
		vg.Length(hmd.Margin)*vg.Inch,
		-vg.Length(hmd.Margin)*vg.Inch)

	title := plot.New()
	title.Title.Text = hmd.Title
	title.Title.TextStyle.Font.Size = vg.Length(hmd.TitleSize) * vg.Inch
	title.HideAxes()
	title.Draw(drc)

	// Compute offset to place the legend
	r := leg.Rectangle(drc)
	legendWidth := r.Max.X - r.Min.X
	legendHeight := r.Max.Y - r.Min.Y

	// Center the legend
	leg.YOffs = (drc.Max.Y - drc.Min.Y - legendHeight) / 2.0

	// Plot the legend and crop the figure
	leg.Draw(drc)
	drc = draw.Crop(drc, 0, -1.5*legendWidth, 0, -1.5*legendWidth)
	plo.Draw(drc)

	// Write out
	f, err := os.Create(hmd.Output)
	if err != nil {
		return err
	}
	defer f.Close()

	_, err = img.WriteTo(f)
	if err != nil {
		return err
	}

	err = f.Close()
	if err != nil {
		return err
	}

	return nil
}
