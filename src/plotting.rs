use crate::params::Derived;
use anyhow::Result;
use ndarray::Array1;
use plotters::coord::Shift;
use plotters::prelude::*;
use std::f64::consts::PI;
use std::sync::Once;

static FONT_INIT: Once = Once::new();

fn init_font() {
    FONT_INIT.call_once(|| {
        let candidates = [
            "C:/Windows/Fonts/arial.ttf",
            "C:/Windows/Fonts/ARIAL.TTF",
            "C:/Windows/Fonts/segoeui.ttf",
            "C:/Windows/Fonts/calibri.ttf",
        ];
        for path in candidates.iter() {
            if let Ok(data) = std::fs::read(path) {
                let leaked: &'static [u8] = Box::leak(data.into_boxed_slice());
                let _ = plotters::style::register_font(
                    "sans-serif",
                    plotters::style::FontStyle::Normal,
                    leaked,
                );
                break;
            }
        }
    });
}

const BLUE: RGBColor = RGBColor(26, 111, 175);
const RED: RGBColor = RGBColor(192, 57, 43);
const GRID_GRAY: RGBColor = RGBColor(200, 200, 200);
const ORANGE: RGBColor = RGBColor(240, 140, 20);

fn stats(x: &Array1<f64>) -> (f64, f64, f64, f64) {
    let n = x.len() as f64;
    let mean = x.iter().sum::<f64>() / n;
    let var = x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
    let std = var.sqrt();
    let max = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min = x.iter().cloned().fold(f64::INFINITY, f64::min);
    (mean, std, min, max)
}

fn viridis(t: f64) -> RGBColor {
    let t = t.clamp(0.0, 1.0);
    let r = (0.267 + t * (0.993 - 0.267)).clamp(0.0, 1.0);
    let g = (0.005 + t * (0.906 - 0.005)).clamp(0.0, 1.0);
    let b = (0.329 + t * (0.144 - 0.329) + (t - 0.5).abs() * 0.4).clamp(0.0, 1.0);
    RGBColor((r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8)
}

fn draw_time_series<DB: DrawingBackend>(
    area: &DrawingArea<DB, Shift>,
    title: &str,
    x_label: &str,
    y_label: &str,
    x: &Array1<f64>,
    y: &Array1<f64>,
    color: RGBColor,
    show_band: bool,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let (mean, std, ymin_v, ymax_v) = stats(y);
    let xmin = x.iter().cloned().fold(f64::INFINITY, f64::min);
    let xmax = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let pad = (ymax_v - ymin_v).abs().max(1e-6) * 0.1;
    let ymin = ymin_v - pad;
    let ymax = ymax_v + pad;
    let mut chart = ChartBuilder::on(area)
        .caption(title, ("sans-serif", 22).into_font())
        .margin(10)
        .x_label_area_size(45)
        .y_label_area_size(70)
        .build_cartesian_2d(xmin..xmax, ymin..ymax)
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc(y_label)
        .label_style(("sans-serif", 16))
        .axis_desc_style(("sans-serif", 18))
        .light_line_style(GRID_GRAY.mix(0.35))
        .draw()
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    if show_band {
        let upper = mean + std;
        let lower = mean - std;
        chart
            .draw_series(std::iter::once(Rectangle::new(
                [(xmin, lower), (xmax, upper)],
                color.mix(0.12).filled(),
            )))
            .map_err(|e| anyhow::anyhow!("{}", e))?;
    }
    chart
        .draw_series(LineSeries::new(
            x.iter().zip(y.iter()).map(|(&a, &b)| (a, b)),
            color.stroke_width(1),
        ))
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    chart
        .draw_series(LineSeries::new(
            [(xmin, mean), (xmax, mean)].into_iter(),
            BLACK.mix(0.6).stroke_width(1),
        ))
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    Ok(())
}

fn draw_phase<DB: DrawingBackend>(
    area: &DrawingArea<DB, Shift>,
    h: &Array1<f64>,
    p_deg: &Array1<f64>,
    t_min: &Array1<f64>,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let hmin = h.iter().cloned().fold(f64::INFINITY, f64::min);
    let hmax = h.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let pmin = p_deg.iter().cloned().fold(f64::INFINITY, f64::min);
    let pmax = p_deg.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let tmin = t_min[0];
    let tmax = t_min[t_min.len() - 1];
    let hpad = (hmax - hmin).abs().max(1e-6) * 0.1;
    let ppad = (pmax - pmin).abs().max(1e-9) * 0.1;
    let mut chart = ChartBuilder::on(area)
        .caption("Faz Portresi (Heave-Pitch)", ("sans-serif", 22).into_font())
        .margin(10)
        .x_label_area_size(45)
        .y_label_area_size(70)
        .build_cartesian_2d(hmin - hpad..hmax + hpad, pmin - ppad..pmax + ppad)
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    chart
        .configure_mesh()
        .x_desc("Heave xi3 (m)")
        .y_desc("Pitch xi5 (derece)")
        .label_style(("sans-serif", 16))
        .axis_desc_style(("sans-serif", 18))
        .light_line_style(GRID_GRAY.mix(0.35))
        .draw()
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    let range = (tmax - tmin).max(1e-9);
    let mut points = Vec::new();
    let mut i = 0usize;
    while i < h.len() {
        let t_norm = ((t_min[i] - tmin) / range).clamp(0.0, 1.0);
        points.push((h[i], p_deg[i], viridis(t_norm)));
        i += 4;
    }
    chart
        .draw_series(
            points
                .into_iter()
                .map(|(x, y, c)| Circle::new((x, y), 1, c.filled())),
        )
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    Ok(())
}

fn draw_psd<DB: DrawingBackend>(
    area: &DrawingArea<DB, Shift>,
    title: &str,
    y_label: &str,
    om: &Array1<f64>,
    s: &Array1<f64>,
    color: RGBColor,
    wn: f64,
    wp: f64,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let mut pairs: Vec<(f64, f64)> = om
        .iter()
        .zip(s.iter())
        .filter(|(&w, &v)| w > 0.0 && w <= 1.0 && v > 0.0)
        .map(|(&w, &v)| (w, v))
        .collect();
    if pairs.is_empty() {
        pairs.push((1e-3, 1e-12));
    }
    let ymin = pairs.iter().map(|p| p.1).fold(f64::INFINITY, f64::min) * 0.5;
    let ymax = pairs.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max) * 2.0;
    let ymin = ymin.max(1e-14);
    let mut chart = ChartBuilder::on(area)
        .caption(title, ("sans-serif", 22).into_font())
        .margin(10)
        .x_label_area_size(45)
        .y_label_area_size(80)
        .build_cartesian_2d(0.0_f64..1.0_f64, (ymin..ymax).log_scale())
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    chart
        .configure_mesh()
        .x_desc("omega (rad/s)")
        .y_desc(y_label)
        .label_style(("sans-serif", 16))
        .axis_desc_style(("sans-serif", 18))
        .light_line_style(GRID_GRAY.mix(0.35))
        .draw()
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    chart
        .draw_series(LineSeries::new(pairs.into_iter(), color.stroke_width(2)))
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    chart
        .draw_series(LineSeries::new(
            [(wn, ymin), (wn, ymax)].into_iter(),
            GREEN.stroke_width(2),
        ))
        .map_err(|e| anyhow::anyhow!("{}", e))?
        .label(format!("wn = {:.3}", wn))
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], GREEN.stroke_width(2)));
    let orange = ORANGE;
    chart
        .draw_series(LineSeries::new(
            [(wp, ymin), (wp, ymax)].into_iter(),
            orange.stroke_width(2),
        ))
        .map_err(|e| anyhow::anyhow!("{}", e))?
        .label(format!("wp = {:.3}", wp))
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], ORANGE.stroke_width(2)));
    chart
        .configure_series_labels()
        .label_font(("sans-serif", 14))
        .background_style(WHITE.mix(0.85))
        .border_style(BLACK)
        .draw()
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    Ok(())
}

pub fn generate_figure(
    out_path: &str,
    t_min: &Array1<f64>,
    heave: &Array1<f64>,
    pitch_rad: &Array1<f64>,
    om3: &Array1<f64>,
    s3: &Array1<f64>,
    om5: &Array1<f64>,
    s5: &Array1<f64>,
    d: &Derived,
) -> Result<()> {
    init_font();
    let root = BitMapBackend::new(out_path, (1950, 1500)).into_drawing_area();
    root.fill(&WHITE).map_err(|e| anyhow::anyhow!("{}", e))?;
    let root = root
        .titled(
            "Eslesik Heave-Pitch Hareketi - GoM Ekstrem Kosulu | Hs = 12.2 m, Tp = 14 s | Kafes Spar Platformu",
            ("sans-serif", 30).into_font(),
        )
        .map_err(|e| anyhow::anyhow!("{}", e))?;
    let rows = root.split_evenly((3, 1));
    let top = &rows[0];
    let mid = rows[1].split_evenly((1, 2));
    let bot = rows[2].split_evenly((1, 2));

    draw_time_series(
        top,
        "Heave Zaman Serisi",
        "Zaman (dk)",
        "Heave xi3 (m)",
        t_min,
        heave,
        BLUE,
        true,
    )?;
    let pitch_deg = pitch_rad.mapv(|v| v * 180.0 / PI);
    draw_time_series(
        &mid[0],
        "Pitch Zaman Serisi",
        "Zaman (dk)",
        "Pitch xi5 (derece)",
        t_min,
        &pitch_deg,
        RED,
        false,
    )?;
    draw_phase(&mid[1], heave, &pitch_deg, t_min)?;
    draw_psd(
        &bot[0],
        "Heave Guc Spektral Yogunlugu",
        "PSD (m^2.s/rad)",
        om3,
        s3,
        BLUE,
        d.wn3,
        d.wp,
    )?;
    draw_psd(
        &bot[1],
        "Pitch Guc Spektral Yogunlugu",
        "PSD (rad^2.s/rad)",
        om5,
        s5,
        RED,
        d.wn5,
        d.wp,
    )?;
    root.present().map_err(|e| anyhow::anyhow!("{}", e))?;
    Ok(())
}
