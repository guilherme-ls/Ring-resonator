import gdsfactory as gf
from gdsfactory.generic_tech import get_generic_pdk, LAYER
from gdsfactory.technology import LayerLevel, LayerStack, LogicalLayer

gf.config.rich_output()
generic_pdk = get_generic_pdk()

# Layer thickness
t_si = 0.220
t_box = 2.0
t_clad = 2.0
t_substrate = 2.0

# Choose layer material and thickness
layer_stack = LayerStack(
    layers={
        "substrate": LayerLevel(
            layer=LogicalLayer(layer=LAYER.WAFER),
            thickness=t_substrate,
            zmin=-t_box-t_substrate,
            material="si",
            background=True,
        ),

        "box": LayerLevel(
            layer=LogicalLayer(layer=LAYER.WAFER),
            thickness=t_box,
            zmin=-t_box,
            material="sio2",
            background=True,
        ),

        "core": LayerLevel(
            layer=LogicalLayer(layer=LAYER.WG),
            thickness=t_si,
            zmin=0.0,
            material="si",
            sidewall_angle=0,
        ),

        "cladding": LayerLevel(
            layer=LogicalLayer(layer=LAYER.WAFER),
            thickness=t_clad,
            zmin=t_si,
            material="sio2",
            background=True,
        ),
    }
)

# Choose waveguide profile
def soi_strip(width=0.450, radius=10.0):
    return gf.cross_section.strip(
        width=width,
        layer=LAYER.WG,
        radius=radius,
    )

# Create PDK
PDK = gf.Pdk(
    name="SOI_220nm",
    base_pdks=[generic_pdk],
    layers=LAYER,
    cross_sections={
        "soi_strip": soi_strip,
    },
    layer_stack=layer_stack,
    layer_views=generic_pdk.layer_views,
)