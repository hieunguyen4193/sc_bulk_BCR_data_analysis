---
title: "Data analysis first dataset - GEX, 26.10.2023"
author:
  - "trnguyen@ukaachen.de"
date: "Last update on 30 October, 2024"
output:
  html_document:
    df_print: paged
    keep_md: true
    toc: true
    toc_float:
      toc_collapsed: false
    toc_depth: 3
    number_sections: false
    theme: lumen
---

<style type="text/css">
script src = "https://ajax.googleapis.com/ajax/libs/jquery/3.4.1/jquery.min.js"
</style>

<script>
 $(document).ready(function() {
    $('body').prepend('<div class=\"zoomDiv\"><img src=\"\" class=\"zoomImg\"></div>');
    // onClick function for all plots (img's)
    $('img:not(.zoomImg)').click(function() {
      $('.zoomImg').attr('src', $(this).attr('src')).css({width: '100%'});
      $('.zoomDiv').css({opacity: '1', width: 'auto', border: '1px solid white', borderRadius: '5px', position: 'fixed', top: '50%', left: '50%', marginRight: '-50%', transform: 'translate(-50%, -50%)', boxShadow: '0px 0px 50px #888888', zIndex: '50', overflow: 'auto', maxHeight: '100%'});
    });
    // onClick function for zoomImg
    $('img.zoomImg').click(function() {
      $('.zoomDiv').css({opacity: '0', width: '0%'}); 
    });
  });
</script>

<style type="text/css">
    div.datatables { height: auto !important;}
</style>




# UMAP

## All samples, grouped by clusters
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

## All samples, highlighted cells with clone VDJ information
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

# VDJ information


## Distribution of mutation rates in all samples

![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

## Mutation rate on UMAP
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

## Violin plot mutation rate
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

## Clone size vs number of mutations
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

## Clone size (>= 10) vs number of mutations
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-9-1.png)<!-- -->



# Differential gene expression analysis between low vs. high mutation rate clones/cells


## UMAP: high mutation rate cells
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-12-1.png)<!-- -->


![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

## DGE table

```{=html}
<div id="htmlwidget-f843a680a1fcd62a0dd1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f843a680a1fcd62a0dd1">{"x":{"filter":"top","vertical":false,"filterHTML":"<tr>\n  <td><\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\" disabled=\"\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0\" data-max=\"1e-15\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0.594310376155458\" data-max=\"0.594310376155459\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\" disabled=\"\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0\" data-max=\"1\" data-scale=\"3\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\" disabled=\"\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0\" data-max=\"1\" data-scale=\"3\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"3.91e-13\" data-max=\"3.92e-13\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","extensions":["Buttons"],"data":[["1"],["Gm30211"],[3.08123707469945e-17],[0.594310376155459],[0.777],[0.566],[3.91163046633095e-13]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Gene<\/th>\n      <th>p_val<\/th>\n      <th>avg_log2FC<\/th>\n      <th>pct.1<\/th>\n      <th>pct.2<\/th>\n      <th>p_val_adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Blfrtip","buttons":["copy","csv","excel","pdf","print"],"lengthMenu":[[10,25,50,-1],["10","25","50","All"]],"columnDefs":[{"targets":"_all","render":"function(data, type, row, meta) {\nreturn type === 'display' && data != null && data.length > 100 ?\n'<span title=\"' + data + '\">' + data.substr(0, 100) + '...<\/span>' : data;\n}"},{"className":"dt-right","targets":[2,3,4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true}},"evals":["options.columnDefs.0.render"],"jsHooks":[]}</script>
```

# Differential gene expression analysis between YFP+ and YFP- cells


```{=html}
<div id="htmlwidget-ba32e0541469b5ed4906" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-ba32e0541469b5ed4906">{"x":{"filter":"top","vertical":false,"filterHTML":"<tr>\n  <td><\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"5.9482839e-08\" data-max=\"3.65741564e-06\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"-0.764113434917828\" data-max=\"0.29101127941599\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0.781\" data-max=\"1\" data-scale=\"3\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0.601\" data-max=\"1\" data-scale=\"3\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"0.000755134644075\" data-max=\"0.046430891547026\" data-scale=\"15\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n<\/tr>","extensions":["Buttons"],"data":[["1","2","3","4"],["Gm42418","Vpreb3","Erp29","mt-Atp8"],[5.94828392339975e-08,1.20594503774938e-07,1.17297943510071e-06,3.65741563978142e-06],[-0.764113434917828,0.291011279415989,0.258246761277864,-0.407189803969179],[1,0.781,0.877,1],[1,0.601,0.733,0.99],[0.000755134644075599,0.00153094722542284,0.0148909739286036,0.0464308915470251]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Gene<\/th>\n      <th>p_val<\/th>\n      <th>avg_log2FC<\/th>\n      <th>pct.1<\/th>\n      <th>pct.2<\/th>\n      <th>p_val_adj<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"Blfrtip","buttons":["copy","csv","excel","pdf","print"],"lengthMenu":[[10,25,50,-1],["10","25","50","All"]],"columnDefs":[{"targets":"_all","render":"function(data, type, row, meta) {\nreturn type === 'display' && data != null && data.length > 100 ?\n'<span title=\"' + data + '\">' + data.substr(0, 100) + '...<\/span>' : data;\n}"},{"className":"dt-right","targets":[2,3,4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true}},"evals":["options.columnDefs.0.render"],"jsHooks":[]}</script>
```

## Highlight YFP cells
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

# Check expression of some genes



## UMAP {.tabset}
### Gene Cxcr3 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-18-1.png)<!-- -->
 
### Gene Tigit 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-18-2.png)<!-- -->
 
### Gene Ly6a 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-18-3.png)<!-- -->
 
### Gene Cd93 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-18-4.png)<!-- -->
 
### Gene H2-Aa 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-18-5.png)<!-- -->
 
### Gene Epcam 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-18-6.png)<!-- -->
 
### Gene Slamf6 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-18-7.png)<!-- -->
 
### Gene Sdc1 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-18-8.png)<!-- -->
 

## Violin plot {.tabset}
### Gene Cxcr3 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-19-1.png)<!-- -->
 
### Gene Tigit 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-19-2.png)<!-- -->
 
### Gene Ly6a 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-19-3.png)<!-- -->
 
### Gene Cd93 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-19-4.png)<!-- -->
 
### Gene H2-Aa 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-19-5.png)<!-- -->
 
### Gene Epcam 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-19-6.png)<!-- -->
 
### Gene Slamf6 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-19-7.png)<!-- -->
 
### Gene Sdc1 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-19-8.png)<!-- -->
 


## Violin plot, YFP vs non YFP {.tabset}
### Gene Cxcr3 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-20-1.png)<!-- -->
 
### Gene Tigit 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-20-2.png)<!-- -->
 
### Gene Ly6a 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-20-3.png)<!-- -->
 
### Gene Cd93 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-20-4.png)<!-- -->
 
### Gene H2-Aa 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-20-5.png)<!-- -->
 
### Gene Epcam 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-20-6.png)<!-- -->
 
### Gene Slamf6 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-20-7.png)<!-- -->
 
### Gene Sdc1 
![](NEW02_DGE_YFP_cells_2nd_dataset_files/figure-html/unnamed-chunk-20-8.png)<!-- -->
 
