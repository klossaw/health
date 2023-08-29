"""
    toolkits for plotting
"""
import warnings
from tqdm import tqdm
from scipy import stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go

from matplotlib import font_manager
from utils import parse_dict_with_default

warnings.filterwarnings("ignore")

PY_PATH = "/cluster/home/bqhu_jh/share/miniconda3/envs/cuda1.7/lib/python3.8/site-packages"
MPL_TTY_PATH = "matplotlib/mpl-data/fonts/ttf"
font_files = font_manager.findSystemFonts(
    fontpaths=f"{PY_PATH}/{MPL_TTY_PATH}", fontext='ttf')
font_manager.fontManager.addfont(f'{PY_PATH}/{MPL_TTY_PATH}/Arial.ttf')
font_manager.fontManager.addfont(f'{PY_PATH}/{MPL_TTY_PATH}/Arial Bold.ttf')

plt.rcParams['font.family'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("white", {"font.family": ["Arial"]})
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# defaults
scatter_x_default = "month"
scatter_order_default = [11, 12, 1, 2, 3, 4, 5, 6]
scatter_hue_default = "period"
scatter_hue_order_default = ["Control-2021", "Control-2022", "Test-2023"]
scatterplus_facet_hue_x_default = "gender"
scatterplus_facet_hue_x_order_default = ["female", "male"]
scatterplus_facet_hue_y_default = "age_groups"
scatterplus_facet_hue_y_order_default = ["<30", "30-45", "45-60", ">60"]

my_pal = ["#0172B6", "#E18727", "#BD3C29", "#21854F", "#7876B1", "#6F99AD", "#00A087", "#EE4C97"]
my_linestyle = [(0, (3, 1, 1, 1)), (0, (5, 10)), "solid", ]
dict_month = {
    1: "Jan", 2: "Feb", 3: "Mar", 4: "Apr",
    5: "May", 6: "Jun", 7: "July", 8: "Aug",
    9: "Sep", 10: "Oct", 11: "Nov", 12: "Dec",
}


def func_q5(value):
    """
        parse np.percentile 5%
    """
    return np.percentile(value, 0.05)


def func_q95(value):
    """
        parse np.percentile 95%
    """
    return np.percentile(value, 0.95)


def get_subplots_with_flank_ratio(n_pos, flank_ratio=0.1):
    """get subplots' coords with flank_ratio

    Args:
        n_pos (int): number of positions
        flank_ratio (float, optional): flank region ratio. Defaults to 0.1.

    Returns:
        np.array([[beg1, end1], [beg2, end2]]): _description_
    """
    n_flanks = n_pos + 1
    flank_len = flank_ratio / n_flanks
    draw_len = (1-flank_ratio) / n_pos
    l_out = []
    for i in range(n_pos):
        beg = flank_len + (draw_len+flank_len)*i
        end = beg + draw_len

        l_out.append([beg, end])

    return np.array(l_out)[::-1]


def get_x_pos_with_n_hues(n_pos, n_hues):
    """get hues' x-coords with number of x and n_hue for each x

    Args:
        n_pos (int): number of positions
        n_hues (float): number of hues for each pos

    Returns:
        np.array([
            [pos1_hue_1, pos1_hue_2, ... , pos1_hue_n],
            [pos2_hue_1, pos2_hue_2, ... , pos2_hue_n]
        ])
    """
    l_x_out = []
    for i in range(n_pos):
        for x_delta in list(np.linspace(-0.5, 0.5, n_hues+2)[1:-1]):
            p_center = i + x_delta
            l_x_out.append(p_center)

    return np.array(l_x_out)


def _generate_x_pos(n_pos, n_hues, hue_ratio=0.15):
    l_x_out = []
    for i in range(n_pos):
        for x_delta in list(np.linspace(-0.5, 0.5, n_hues+2)[1:-1]):
            p_delta = 1/(n_hues+2)
            p_center = i + 1 + x_delta
            l_x_out.append(p_center-hue_ratio*p_delta)
            l_x_out.append(p_center+hue_ratio*p_delta)

    return np.array(l_x_out)

def _parse_column_state(item, l_tags, default_result="Normal"):
    for tag in l_tags:
        if item[tag] > 0:
            return tag

    return default_result

def arr_ratio_pval(arr1, arr2):
    """calculate t-test [fold-change and pvalue]

    Args:
        arr1 (np.array(np.float))): array1 to compare
        arr2 (np.array(np.float))): array2 to compare

    Returns:
        [float, float]: fold-change and -log10(pvalue)
    """
    ratio_t = np.mean(arr1) / (np.mean(arr2)+1e-5)
    if ratio_t > 2:
        ratio_t = 2

    if ratio_t < 0.5:
        ratio_t = 0.5

    _, p_value = stats.ttest_ind(arr1, arr2)
    p_value_log10 = -1*np.log10(p_value)
    return [ratio_t, p_value_log10]

class Figure(object):
    """plot with Figure objs

    Including:
        Scatter
        Stack
        Sankey
    """
    def __init__(self, figsize=(6,6), rename_dict=None, n_cols=4, n_rows=4):
        self.figsize = figsize
        self.rename_dict = rename_dict
        self.n_rows = n_rows
        self.n_cols = n_cols


class Scatter(Figure):
    """My scatter plot for df_table1plus

    Supports:
        plot_scatter_ax
        plot_scatter
        plot_scatter_ax_plus
    """
    def __init__(self, figsize=(6, 6), rename_dict=None, n_cols=4, n_rows=4):
        super().__init__(figsize, rename_dict, n_cols, n_rows)

    def plot_scatter_ax(self, df, y, axes,
                        x=scatter_hue_default, order=None, hue=scatter_hue_default, 
                        hue_order=None, show_error_bar=False):
        """plot scatter plot for one tags by month, given ax.

        Args:
            df (_type_): _description_
            y (_type_): _description_
            axes (_type_): _description_
            x (_type_, optional): _description_. Defaults to scatter_hue_default.
            order (_type_, optional): _description_. Defaults to None.
            hue (_type_, optional): _description_. Defaults to scatter_hue_default.
            hue_order (_type_, optional): _description_. Defaults to None.
            show_error_bar (bool, optional): _description_. Defaults to False.

        Returns:
            _type_: _description_
        """
        if order is None:
            order = list(df[x].drop_duplicates())
            if scatter_order_default is not None:
                order = scatter_order_default

        if hue_order is None:
            hue_order = list(df[hue].drop_duplicates())
            if scatter_hue_order_default is not None:
                hue_order = scatter_hue_order_default

        df_my_pvt = df[[x, y, hue]].dropna().pivot_table(
            index=hue, columns=x, values=y,
            aggfunc=[len, np.mean, np.std, func_q5, func_q95]
        )

        df_mean_sub = df_my_pvt["mean"][order]
        df_n_sub = df_my_pvt["len"][order]
        df_q5_sub = df_my_pvt["func_q5"][order]
        df_q95_sub = df_my_pvt["func_q95"][order]

        n_points = len(order)
        pos_move = np.linspace(-0.2, 0.2, len(hue_order))
        for idx,label in enumerate(hue_order):
            axes.plot(np.arange(n_points)+pos_move[idx],
                        df_mean_sub.loc[label],
                        linestyle=my_linestyle[idx],
                        color=my_pal[idx], label=label
            )
            scatter_obj = axes.scatter(np.arange(n_points)+pos_move[idx],
                        df_mean_sub.loc[label],
                        s=df_n_sub.loc[label]/100, color=my_pal[idx]
            )
            if show_error_bar:
                axes.errorbar(np.arange(n_points)+pos_move[idx],
                        df_mean_sub.loc[label],
                        yerr=[df_q5_sub.loc[label], df_q95_sub.loc[label]],
                        color=my_pal[idx]
                )

        axes.axvspan(0.5, 2.5, alpha=0.1, color='blue')
        name = parse_dict_with_default(y, self.rename_dict)
        axes.set_title(f"{name}")
        df_my_pvt["tag"] = y
        return scatter_obj, df_my_pvt


    def plot_scatter(self, df, cols,
                     x=scatter_x_default, order=None,
                     hue=scatter_hue_default, hue_order=None, show_error_bar=False,
                    ylim=None):
        """scatter plots for df hue by cols

        Args:
            df (_type_): _description_
            cols (_type_): _description_
            x (_type_, optional): _description_. Defaults to scatter_x_default.
            order (_type_, optional): _description_. Defaults to None.
            hue (_type_, optional): _description_. Defaults to scatter_hue_default.
            hue_order (_type_, optional): _description_. Defaults to None.
            show_error_bar (bool, optional): _description_. Defaults to False.
            ylim (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """
        if order is None:
            order = list(df[x].drop_duplicates())
            if scatter_order_default is not None:
                order = scatter_order_default

        if hue_order is None:
            hue_order = list(df[hue].drop_duplicates())
            if scatter_hue_order_default is not None:
                hue_order = scatter_hue_order_default


        l_dfs = []
        fig = plt.figure(figsize=self.figsize)
        for i in tqdm(range(len(cols))):
            if i >= self.n_cols*self.n_rows:
                break

            axes = fig.add_subplot(self.n_rows, self.n_cols, 1+i)
            tag = cols[i]
            scatter_obj, df_my_pvt = self.plot_scatter_ax(df, x=x,
                                    y=tag, axes=axes, order=order, hue=hue,
                                    hue_order=hue_order, show_error_bar=show_error_bar)
            l_dfs.append(df_my_pvt)

            if i > len(cols) - self.n_rows - 1:
                axes.set_xticks(np.arange(len(order)))
                axes.set_xticklabels([parse_dict_with_default(mon, dict_month)
                                for mon in order], rotation=305)
            else:
                axes.set_xticks([])
                axes.set_xticklabels([])

            if ylim is not None:
                axes.set_ylim(ylim)

            if i % self.n_cols == self.n_cols-1:
                if i < self.n_cols:
                    axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                    continue

                kwargs = dict(prop="sizes", num=5, fmt="{x:.0f}", func=lambda x: x*50)
                axes.legend(*scatter_obj.legend_elements(**kwargs),
                        loc='center left', bbox_to_anchor=(1, 0.5))

        return fig, pd.concat(l_dfs)


    def plot_scatter_ax_plus(self, df, y, x=scatter_x_default, order=None,
                             hue=scatter_hue_default, hue_order=None,
                             facet_hue_x="gender", facet_hue_x_order=None,
                             facet_hue_y="age_groups", facet_hue_y_order=None,
                             show_error_bar=False):
        """plot scatter plot by age-gender

        Args:
            df (_type_): _description_
            y (_type_): _description_
            x (_type_, optional): _description_. Defaults to scatter_x_default.
            order (_type_, optional): _description_. Defaults to None.
            hue (_type_, optional): _description_. Defaults to scatter_hue_default.
            hue_order (_type_, optional): _description_. Defaults to None.
            facet_hue_x (str, optional): _description_. Defaults to "gender".
            facet_hue_x_order (_type_, optional): _description_. Defaults to None.
            facet_hue_y (str, optional): _description_. Defaults to "age_groups".
            facet_hue_y_order (_type_, optional): _description_. Defaults to None.
            show_error_bar (bool, optional): _description_. Defaults to False.

        Returns:
            _type_: _description_
        """
        if order is None:
            order = list(df[x].drop_duplicates())
            if scatter_order_default is not None:
                order = scatter_order_default

        if hue_order is None:
            hue_order = list(df[hue].drop_duplicates())
            if scatter_hue_order_default is not None:
                hue_order = scatter_hue_order_default

        if facet_hue_x_order is None:
            facet_hue_x_order = list(df[facet_hue_x].drop_duplicates())
            if scatterplus_facet_hue_x_order_default is not None:
                facet_hue_x_order = scatterplus_facet_hue_x_order_default

        if facet_hue_y_order is None:
            facet_hue_y_order = list(df[facet_hue_x].drop_duplicates())
            if scatterplus_facet_hue_y_order_default is not None:
                facet_hue_y_order = scatterplus_facet_hue_y_order_default

        df_my_pvt = df[[y, facet_hue_x, facet_hue_y, hue, x]].dropna().pivot_table(
            index=[facet_hue_x, facet_hue_y, hue], columns=x,
            aggfunc=[len, np.mean, np.std, func_q5, func_q95], values=y
        )
        np_val = df_my_pvt["mean"].values
        np_min, np_max = np_val.min(), np_val.max()
        fig = plt.figure(figsize=self.figsize)
        for idx_i, gender in enumerate(facet_hue_x_order):
            for idx, age_group in enumerate(facet_hue_y_order):
                axes = fig.add_subplot(len(facet_hue_x_order),
                                       len(facet_hue_y_order), idx_i*len(facet_hue_y_order)+idx+1)
                df_table1plus_sub = df[
                    (df[facet_hue_x] == gender) &
                    (df[facet_hue_y] == age_group)
                ]

                scatter_obj, _ = self.plot_scatter_ax(df_table1plus_sub, x=x,
                                            y=y, axes=axes, order=order, hue=hue,
                                            hue_order=hue_order, show_error_bar=show_error_bar)

                axes.set_title(f"{gender}_{age_group}")
                axes.set_ylim([np_min, np_max])
                if idx > 0:
                    axes.set_yticks([])
                    axes.set_yticklabels([])

                if idx_i == 0:
                    axes.set_xticks([])
                    axes.set_xticklabels([])
                    if idx == 0:
                        name = parse_dict_with_default(y, self.rename_dict)
                        axes.set_ylabel(f"{name}")
                    if idx == 3:
                        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))

                if idx_i == 1:
                    axes.set_xticks(np.arange(len(order)))
                    axes.set_xticklabels([parse_dict_with_default(
                        mon, dict_month) for mon in order], rotation=305)
                    if idx == 3:
                        kwargs = dict(prop="sizes", num=5,
                                fmt="{x:.0f}", func=lambda x: x*100)
                        axes.legend(*scatter_obj.legend_elements(**kwargs),
                                loc='center left', bbox_to_anchor=(1, 0.5))

        return fig, df_my_pvt


class Stack(Figure):
    """My stack plot for df_table1plus

    Supports:
        plot_bar_stacked
    """
    def __init__(self, figsize=(6, 6), rename_dict=None, n_cols=4, n_rows=4):
        super().__init__(figsize, rename_dict, n_cols, n_rows)

    def _pre_process(self, df_data, x_value, hue, y_value, y_order):
        df_tmp = df_data[[x_value, hue]]
        df_tmp[y_value] = df_data.apply(
            lambda line: _parse_column_state(line, y_order), axis=1)

        df_tmp_grouped = df_tmp.groupby([hue, x_value, y_value]).size().unstack().reset_index()

        # 转换为百分比
        df_tmp_grouped['total'] = df_tmp_grouped.sum(axis=1)
        df_tmp_grouped[y_order + ["Normal"]] = 100*df_tmp_grouped[y_order + ["Normal"]].values / \
            df_tmp_grouped['total'].values.reshape(-1, 1)

        return df_tmp_grouped

    def _plot_stacked_barplot(self, df_grouped, x_value, hue, y_value, y_order,
                            axes, order, hue_order, my_cmap=None):
        if my_cmap is None:
            my_cmap = my_pal

        l_labels = y_order + ["Normal"]
        l_colors = my_cmap
        my_width = 1 / (len(hue_order)+1)
        x_pos_all = get_x_pos_with_n_hues(len(order), len(hue_order)).reshape(len(order), -1)
        for idx, hue_name in enumerate(hue_order):
            x_pos = x_pos_all[:, idx]
            df_grouped_sub = df_grouped[df_grouped[hue] == hue_name]
            y_total = 0
            for idx_label,label in enumerate(l_labels):
                y_plot = np.array([
                    df_grouped_sub[df_grouped_sub[x_value] == x_rank][label].values[0]
                        for x_rank in order
                ])
                kwargs = {"width": my_width,
                        "color": l_colors[idx_label], "bottom": y_total}
                if idx == 0:
                    kwargs["label"] = label

                axes.bar(x_pos, y_plot, **kwargs)
                y_total += y_plot

        axes.set_title(y_value)
        axes.set_xlabel(x_value)
        axes.set_ylabel("Percentage")
        axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        axes.set_xticks(range(len(order)))
        axes.set_xticklabels([parse_dict_with_default(mon, dict_month) for mon in order],
                            ha='center', va="center")


    def plot_bar_stacked(self, df, x, hue, y, y_order, order=None, hue_order=None,
                         cmap=None):
        """plot bar stacked with x,y and hue

        Args:
            df (_type_): _description_
            x (_type_): _description_
            hue (_type_): _description_
            y (_type_): _description_
            y_order (_type_): _description_
            order (_type_, optional): _description_. Defaults to None.
            hue_order (_type_, optional): _description_. Defaults to None.
            cmap (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """
        if cmap is None:
            cmap = my_pal

        if order is None:
            order = sorted(set(df[x]))
        if hue_order is None:
            hue_order = sorted(set(df[hue]))

        fig = plt.figure(figsize=self.figsize)
        axes = fig.add_subplot(1, 1, 1)
        df_grouped = self._pre_process(df, x, hue, y, y_order)
        self._plot_stacked_barplot(df_grouped, x, hue, y, y_order,
                            axes, order, hue_order, cmap)
        return fig, df_grouped


class Sankey(Figure):
    """My sankey plot for df_table1plus

    Args:
        Figure (_type_): _description_
    """
    def __init__(self, figsize=(6, 6), rename_dict=None, n_cols=4, n_rows=4):
        super().__init__(figsize, rename_dict, n_cols, n_rows)

    def _plotly_sankey(self, nodes, edges, domain=None):
        # 创建sankey图
        sankey_plot = go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color='black', width=0.5),
                label=[node['label'] for node in nodes],
            ),
            link=dict(
                source=[edge['source'] for edge in edges],
                target=[edge['target'] for edge in edges],
                value=[edge['value'] for edge in edges],
            ),
            domain=domain
        )

        # 显示图形
        return sankey_plot

    def _get_pair_value(self, df_cnt_sub_pvt, i, j):
        try:
            val = df_cnt_sub_pvt.loc[i, j].values[0]
        except:
            val = 0

        return val

    def _get_pair_edge(self, df_cnt_sub_pvt, idx_pair, i, j):
        val = self._get_pair_value(df_cnt_sub_pvt, i, j)
        return {'source': 2*idx_pair+i, 'target': 2*(idx_pair+1)+j, 'value': val}


    def plot_sankey(self, df_cnt_sub, l_pairs, domain=None):
        """plot sankey from df_cnt_sub

        Args:
            df_cnt_sub (_type_): _description_
            l_pairs (_type_): _description_
            domain (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_


        df_cnt_sub = 
            |      | variable       |   month | age_groups   | gender   |   2021 |   2022 |   2023 |   0 |
            |-----:|:---------------|--------:|:-------------|:---------|-------:|-------:|-------:|----:|
            |   91 | HEART.T_change |       1 | 45-60        | female   |      0 |      0 |      0 | 111 |
            |  452 | HEART.T_change |       1 | 45-60        | female   |      0 |      0 |      1 |  20 |
            |  912 | HEART.T_change |       1 | 45-60        | female   |      0 |      1 |      1 |   7 |
            | 1023 | HEART.T_change |       1 | 45-60        | female   |      1 |      0 |      1 |   5 |
            | 1033 | HEART.T_change |       1 | 45-60        | female   |      1 |      1 |      1 |   5 |
            | 1241 | HEART.T_change |       1 | 45-60        | female   |      0 |      1 |      0 |   3 |
            | 1309 | HEART.T_change |       1 | 45-60        | female   |      1 |      0 |      0 |   2 |
            | 1361 | HEART.T_change |       1 | 45-60        | female   |      1 |      1 |      0 |   1 |
        """
        nodes = [{"label": f"{l_pairs[0][0]}_0"}, {"label": f"{l_pairs[0][0]}_1"}]
        edges = []
        for idx_pair,pair in enumerate(l_pairs):
            nodes.append({"label": f"{pair[1]}_0"})
            nodes.append({"label": f"{pair[1]}_1"})

            df_cnt_sub_pvt = df_cnt_sub[pair+[0]].pivot_table(index=pair, values=0, aggfunc=np.sum)
            edges.append(self._get_pair_edge(df_cnt_sub_pvt, idx_pair, 0, 0))
            edges.append(self._get_pair_edge(df_cnt_sub_pvt, idx_pair, 0, 1))
            edges.append(self._get_pair_edge(df_cnt_sub_pvt, idx_pair, 1, 0))
            edges.append(self._get_pair_edge(df_cnt_sub_pvt, idx_pair, 1, 1))

        return self._plotly_sankey(nodes, edges, domain), {"nodes": nodes, "edges": edges}


    def plot_sankey_subplots(self, df_cnt, month, tag, l_age_groups, l_pairs, prefix="Figure2"):
        """plot multiple sankey plots

        Args:
            df_cnt (_type_): _description_
            month (_type_): _description_
            tag (_type_): _description_
            l_age_groups (_type_): _description_
            l_pairs (_type_): _description_
            prefix (str, optional): _description_. Defaults to "Figure2".

        Returns:
            _type_: _description_
        """
        data = []
        layout =  go.Layout(
            title = f"{tag}, month:{month}",
            font = dict(
            size = 10
            )
        )
        sankey_dicts = {}
        df_cnt_sub = df_cnt[(df_cnt["variable"]==tag) & (df_cnt["month"].isin(month))]
        l_pos = get_subplots_with_flank_ratio(len(l_age_groups), flank_ratio=0.3)
        for idx,age_groups in enumerate(l_age_groups):
            df_cnt_sub_female = df_cnt_sub[(df_cnt_sub["gender"]=="female") &
                                           (df_cnt_sub["age_groups"]==age_groups)]
            domain = {
                    'x': [0, 0.45],
                    'y': l_pos[idx]
            }
            sankey_obj, sankey_dict_female = self.plot_sankey(df_cnt_sub_female, l_pairs, domain)
            data.append(sankey_obj)
            df_cnt_sub_male = df_cnt_sub[(df_cnt_sub["gender"]=="male") &
                                    (df_cnt_sub["age_groups"]==age_groups)]
            domain = {
                    'x': [0.55, 1.0],
                    'y': l_pos[idx]
            }
            sankey_obj, sankey_dict_male = self.plot_sankey(df_cnt_sub_male, l_pairs, domain)
            data.append(sankey_obj)
            sankey_dicts[age_groups] = {"female": sankey_dict_female, "male": sankey_dict_male}

        fig = go.Figure(data=data, layout=layout)
        fig.update_layout(
            autosize=False,
            width=500, height=600,
            margin=dict(
                pad=1
            ),
        )
        return fig, sankey_dicts


class BxxPvalue(Figure):
    """boxplot/barplot with P-value above

    Supports:
        get_fc_pval
        plot_bxxplot_pvalue
    """
    def __init__(self, figsize=(6, 6), rename_dict=None, n_cols=4, n_rows=4):
        super().__init__(figsize, rename_dict, n_cols, n_rows)
        self.hue_order = [
            'Control-2021_<30',   'Control-2022_<30',   'Test-2023_<30', "",
            'Control-2021_30-45', 'Control-2022_30-45', 'Test-2023_30-45', "",
            'Control-2021_45-60', 'Control-2022_45-60', 'Test-2023_45-60', "",
            'Control-2021_>60',   'Control-2022_>60',   'Test-2023_>60',"",
        ]
        self.hue_order_1 = [
            'Control-2021',   'Control-2022',   'Test-2023'
        ]
        self.hue_order_2 = ["<30", "30-45", "45-60", ">60"]

        self.l_months = [11,12, 1,2,3,  4,5,6]
        self.l_genders = ["female", "male"]
        self.l_age_groups = ["<30", "30-45", "45-60", ">60"]


    def get_fc_pval(self, df_table1plus, y, x="month", hue="period_age", hue_col="gender",
                        hue_col_order=None, order=None):
        """get fold-change and pvalue from boxplot

        Args:
            df_table1plus (pd.DataFrame): raw input dataframe
            tag (_type_): feature name

        Returns:
            _type_: _description_
        """
        df_month_var = df_table1plus[[hue, x, hue_col, y]].dropna()
        l_out = []
        for gender in hue_col_order:
            for month in order:
                for age_group in self.hue_order_2:
                    hue_t  = f"Test-2023_{age_group}"
                    hue_c1 = f"Control-2022_{age_group}"
                    hue_c2 = f"Control-2021_{age_group}"
                    df_p_plot = df_month_var[
                                    (df_month_var[hue_col]==gender) &
                                    (df_month_var[x] == month)
                    ]
                    if (df_p_plot is None) or (df_p_plot.shape[0]==0):
                        continue

                    subset_t  = df_p_plot[df_p_plot[hue] == hue_t][y].dropna()
                    subset_c1 = df_p_plot[df_p_plot[hue] == hue_c1][y].dropna()
                    subset_c2 = df_p_plot[df_p_plot[hue] == hue_c2][y].dropna()
                    n_t = min(len(subset_t), len(subset_c1))
                    n_c = min(len(subset_c1), len(subset_c2))
                    fc_t, pv_t = arr_ratio_pval(subset_t, subset_c1)
                    fc_c, pv_c = arr_ratio_pval(subset_c1, subset_c2)

                    l_out.append(
                        pd.Series([y, gender, age_group, month, "2021-2022", n_c, fc_c, pv_c])
                    )
                    l_out.append(
                        pd.Series([y, gender, age_group, month, "2022-2023", n_t, fc_t, pv_t])
                    )

        df_out = pd.DataFrame(l_out)
        df_out.columns = [
            "item_id", "gender", "age_group", "month", 
            "period", "n", "fold_change", "log10_p"
        ]
        return df_out

    def _plot_bxxplot_with_pvalue(self, df_month_var, x, y, order, hue,
                                  hue_order, hue_col, hue_col_order, sns_type, ax1, ax2):
        df_fc_pval = self.get_fc_pval(df_month_var, y, x=x, hue=hue, order=order,
                                                    hue_col=hue_col, hue_col_order=hue_col_order)

        df_tmp = df_fc_pval[(df_fc_pval["item_id"]==y)]
        n_points = len(order)

        sub_idx0 = list(filter(lambda x: x%2==0, range(df_tmp.shape[0])))
        sub_idx1 = list(filter(lambda x: x%2==1, range(df_tmp.shape[0])))
        df_tmp0 = df_tmp.iloc[sub_idx0]
        df_tmp1 = df_tmp.iloc[sub_idx1]

        np_x1 = _generate_x_pos(n_points, len(self.hue_order_1)+1)
        scale_point_size = 50
        minimal_point_size = 2
        max_log10_pvalue = 5
        p_value_0 = df_tmp0["log10_p"].values
        p_value_0[p_value_0>max_log10_pvalue] = max_log10_pvalue
        scatter_obj1 = ax1.scatter(np_x1[sub_idx0], df_tmp0["fold_change"],
                    s=df_tmp0["n"]/scale_point_size+minimal_point_size,
                    cmap="viridis",
                    c=p_value_0,
                    marker="^",
                    vmin=0, vmax=max_log10_pvalue,
                    label="fold-change_22-21"
        )
        p_value_1 = df_tmp1["log10_p"].values
        p_value_1[p_value_1>max_log10_pvalue] = max_log10_pvalue
        np_x2 = _generate_x_pos(n_points, 4)
        scatter_obj2 = ax1.scatter(np_x2[sub_idx1], df_tmp1["fold_change"],
                    s=df_tmp1["n"]/scale_point_size+minimal_point_size,
                    cmap="viridis",
                    c=p_value_1,
                    marker="o",
                    vmin=0, vmax=max_log10_pvalue,
                    label="fold-change_23-22"
        )

        ax1.hlines(y=1, xmin=-1, xmax=n_points+2, colors="gray", linestyles="-", linewidth=0.5)
        ax1.set_xticks(range(1, n_points+1))
        ax1.set_xticklabels(order)
        ax1.set_xticklabels([])
        ax1.set_xlim(0.5, n_points+0.5)
        y_max = max(df_tmp0["fold_change"].max(), df_tmp1["fold_change"].max()) * 1.05
        y_min = min(df_tmp0["fold_change"].min(), df_tmp1["fold_change"].min()) * 0.95
        y_max = max(y_max, 1.1)
        y_min = min(y_min, 0.9)
        ax1.set_ylim(y_min, y_max)

        # shape -> comparision
        legend_shape = ax1.legend(title="comparasion", loc='upper left',
                                    bbox_to_anchor=(-0.34, 0.5))
        ax1.add_artist(legend_shape)

        # color -> pvalue
        sc = scatter_obj1
        if p_value_1.max() > p_value_0.max():
            sc = scatter_obj2

        legend_color = ax1.legend(*sc.legend_elements(num=3), title="log10(p-value)",
                                  loc='center left', bbox_to_anchor=(1.3, 0.5)
        )
        ax1.add_artist(legend_color)

        # size -> n-people
        kwargs = dict(prop="sizes", num=5, fmt="{x:.0f}",
                      func=lambda x: (x-minimal_point_size)*scale_point_size)
        ax1.legend(*scatter_obj2.legend_elements(**kwargs), title="n-people",
                        loc='lower left', bbox_to_anchor=(1, 0.5))

        # if sns_type == "boxplot" or sns_type == "box":
        if sns_type in {"boxplot", "box"}:
            sns.boxplot(df_month_var, x=x, y=y, palette=my_pal[0:4],
                        hue=hue, hue_order=hue_order, order=order,
                        ax=ax2, showfliers = False)

        if sns_type in {"barplot", "bar"}:
            sns.barplot(df_month_var, x=x, y=y, palette=my_pal[0:4],
                        hue=hue, hue_order=hue_order, order=order,
                        ax=ax2)
        ax2.legend_.remove()
        ax2.set_ylabel("")
        return df_fc_pval

    def plot_bxxplot_pvalue(self, df_table1plus, y, x="month", order=None, hue="period_age",
                            hue_order=None, hue_col="gender", hue_col_order=None,
                            sns_type="boxplot"):
        """plot boxplot with pvalue above

        Args:
            df_table1plus (pd.DataFrame): raw input dataframe
            y (_type_): feature name

        Returns:
            _type_: _description_
        """
        if hue_col_order is None:
            hue_col_order = self.l_genders
        if hue_order is None:
            hue_order = self.hue_order
        if order is None:
            order = self.l_months

        fig = plt.figure(figsize=self.figsize)
        ax1 = fig.add_axes([0.15, 0.84, 0.80, 0.15])
        ax2 = fig.add_axes([0.15, 0.55, 0.80, 0.27])

        df_month_var = df_table1plus[["period_age", "month", hue_col, y]].dropna()
        df_fcp1 = self._plot_bxxplot_with_pvalue(
            df_month_var[df_month_var[hue_col]==hue_col_order[0]],
            x, y, order, hue, hue_order, hue_col, hue_col_order, sns_type, ax1, ax2
        )
        ax1.set_title(parse_dict_with_default(y, self.rename_dict))
        ax2.set_xticklabels([])
        ax2.set_xlabel("")
        ax3 = fig.add_axes([0.15, 0.34, 0.80, 0.15])
        ax4 = fig.add_axes([0.15, 0.05, 0.80, 0.27])
        df_fcp2 = self._plot_bxxplot_with_pvalue(
            df_month_var[df_month_var[hue_col]==hue_col_order[1]],
            x, y, order, hue, hue_order, hue_col, hue_col_order, sns_type, ax3, ax4
        )
        df_fc_pval = pd.concat([df_fcp1, df_fcp2])
        return fig, df_fc_pval
